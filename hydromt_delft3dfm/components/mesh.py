"""DFlowFM Mesh Component."""

import glob
import logging
from os.path import dirname, join
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
import xugrid as xu
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import MeshComponent
from pyproj import CRS
from shapely.geometry import box

from hydromt_delft3dfm import workflows
from hydromt_delft3dfm.utils import io_utils, mesh_utils

__all__ = ["DFlowFMMeshComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DFlowFMMeshComponent(MeshComponent):
    """
    Manage the Delft3D FM mesh files for model meshes.

    This class is used to manage unstructured mesh data in a model. The mesh component
    data stored in the ``data`` property is a xugrid.UgridDataset object.
    """

    # Mesh2d resolution
    res: float | None = None

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "dflowfm/fm_net.nc",
    ):
        """Initialize DFlowFMMeshComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            by default "dflowfm/fm_net.nc", i.e. one file per mesh in
            the data dictionary.
        """
        super().__init__(
            model=model,
            filename=filename,
        )

    # TODO: improved hydromt-core method, fix in hydromt-core instead
    @property
    def crs(self) -> CRS | None:
        """Returns model mesh crs."""
        if not self.is_empty:
            return next(iter(self.data.ugrid.crs.values()))
        return None

    # TODO: improved hydromt-core method, fix in hydromt-core instead
    @property
    def bounds(self) -> tuple[float, float, float, float] | None:
        """Returns model mesh bounds."""
        if not self.is_empty:
            return self.data.ugrid.bounds
        return None

    # TODO: improved hydromt-core method, fix in hydromt-core instead
    @property
    def _region_data(self) -> gpd.GeoDataFrame | None:
        """Return mesh total_bounds as a geodataframe."""
        if not self.is_empty:
            region = gpd.GeoDataFrame(
                geometry=[box(*self.data.ugrid.total_bounds)], crs=self.crs
            )
            return region
        return None

    # TODO: add to hydromt-core instead
    @property
    def is_empty(self):
        """Check whether the mesh is empty or not."""
        ndims = len(self.data.sizes)
        if ndims == 0:
            return True
        else:
            return False

    ### I/O methods ###
    @hydromt_step
    def read(self) -> None:
        """Read network file with Hydrolib-core and extract mesh/branches info."""
        self.root.is_reading_mode()

        # Read mesh
        # hydrolib-core convention
        fn_network = self.model.dfmmodel.geometry.netfile.filepath
        if fn_network is None:
            raise ValueError(
                "hydromt_delft3dfm cannot read a model without a mesh/network."
            )
        network = self.model.dfmmodel.geometry.netfile.network

        # read the crs from the network with xugrid. TODO: this should be done by
        # hydrolib-core instead https://github.com/Deltares/HYDROLIB-core/issues/1047
        mdu_folder = self.model.dfmmodel.filepath.parents[0]
        fp_network = join(mdu_folder, fn_network)
        uds = xu.open_dataset(fp_network)
        crs_network = next(iter(uds.ugrid.crs.values()))

        # as a fallback, get it from one of the geoms. Cannot use geoms.read() yet
        # because for some geoms (crosssections, manholes) mesh needs to be read first.
        geoms_fns = glob.glob(join(self.root.path, "geoms", "*.geojson"))
        if len(geoms_fns) > 0:
            crs_geoms = gpd.read_file(geoms_fns[0]).crs
        else:
            crs_geoms = None

        # set the model crs, prefer the CRS from the network, then from the geoms.
        # If both are not found, self.model._crs is not overwritten so the value
        # provided to DFlowFMModel() is used (None per default).
        if crs_network:
            logger.debug("found CRS in the mesh")
            if self.model._crs:
                logger.debug(
                    "CRS is provided, but it is overwritten by the CRS from the "
                    "mesh/network file."
                )
            self.model._crs = crs_network
        elif crs_geoms:
            logger.debug("found CRS in the geoms")
            if self.model._crs:
                logger.debug(
                    "CRS is provided, but it is overwritten by the CRS from the geoms."
                )
            self.model._crs = crs_geoms

        # raise an error if the crs was not found in the mesh, nor in the geoms, nor
        # was provided to DFlowFMModel().
        if not self.model._crs:
            raise ValueError(
                "CRS was not found in the mesh or the geoms of the model, please pass "
                "it as an argument to DFlowFMModel() or as global property in your "
                "workflow yml file."
            )

        crs = self.model.crs

        # convert to xugrid
        mesh = mesh_utils.mesh_from_hydrolib_network(network, crs=crs)
        # set mesh
        self._data = mesh

        # update resolution
        if "mesh2d" in self.mesh_names:
            if self.res is None:
                self.res = np.max(np.diff(self.mesh_grids["mesh2d"].node_x))

        # creates branches geometry from network1d
        if "network1d" in self.mesh_names:
            network1d_dataset = self.mesh_datasets["network1d"]
            # Create the branches GeoDataFrame (from geom)
            # network1d_geometry = self.mesh_gdf["network1d"] this returns the network
            branches = mesh_utils.network1d_geoms_geodataframe(network1d_dataset)
            # branches["branchtype"] = network1d_dataset["network1d_branch_type"]
            # might support in the future
            # https://github.com/Deltares/HYDROLIB-core/issues/561

            # Add branchtype, properties from branches.gui file
            logger.info("Reading branches GUI file")
            branches = io_utils.read_branches_gui(branches, self.model.dfmmodel)

            # Set branches
            self.model.set_branches(branches)

    @hydromt_step
    def write(self, write_gui: bool = True) -> None:
        """Write 1D branches and 2D mesh at <root/dflowfm/fm_net.nc>."""
        self.root.is_writing_mode()
        if self.is_empty:
            raise RuntimeError(
                "hydromt_delft3dfm cannot write a model without a mesh/network."
            )
        logger.info("Writing mesh file.")

        # get mesh savedir from dimr (same as mdu path)
        mdu_filename = self.model.mdu._filename
        fm_workingdir = dirname(mdu_filename)
        savedir = join(self.root.path, fm_workingdir)
        Path(savedir).mkdir(parents=True, exist_ok=True)
        mesh_filename = self.model.mdu.get_value(key="geometry.netfile")
        # fallback argument does not work, so repace None with default value manually
        if mesh_filename is None:
            mesh_filename = "fm_net.nc"

        # write mesh
        # HydroMT convention - FIXME hydromt-core does not seem to read the 1D and
        # links part of the mesh
        # super().write_mesh(fn=join(savedir, mesh_filename))

        # write with hydrolib-core
        # Note: hydrolib-core writes more information including attributes and
        # converts some variables using start_index
        # FIXME: does not write crs that is recongnised by Delft3D FM GUI.
        # check dfm_tools/meshkernel_helpers.py#L82

        network = mesh_utils.hydrolib_network_from_mesh(self.data)
        network.to_file(Path(join(savedir, mesh_filename)))

        # save relative path to mdu
        self.model.mdu.set("geometry.netfile", mesh_filename)

        # other mesh1d related geometry TODO update
        if "mesh1d" in self.mesh_names and write_gui and not self.model.branches.empty:
            logger.info("Writing branches.gui file")
            if "manholes" in self.model.geoms.data:
                io_utils.write_branches_gui(self.model.branches, savedir)

    ### Mutating methods ###
    def set(
        self,
        data: xu.UgridDataArray | xu.UgridDataset,
        name: str | None = None,
        grid_name: str | None = None,
        overwrite_grid: bool = False,
    ) -> None:
        """Add data to mesh.

        All layers of mesh have identical spatial coordinates in Ugrid conventions.

        Parameters
        ----------
        data: xugrid.UgridDataArray or xugrid.UgridDataset
            new layer to add to mesh
        name: str, optional
            Name of new object layer, this is used to overwrite the name of
            a UgridDataArray.
        grid_name: str, optional
            Name of the mesh grid to add data to. If None, inferred from data.
            Can be used for renaming the grid.
        overwrite_grid: bool, optional
            If True, overwrite the grid with the same name as the grid in self.mesh.
        """
        # Check if new grid_name
        if grid_name not in self.mesh_names:
            new_grid = True
        else:
            new_grid = False

        # First call super method to add mesh data
        super().set(
            data=data,
            name=name,
            grid_name=grid_name,
            overwrite_grid=overwrite_grid,
        )

        # check if 1D and 2D and and 1D2D links and overwrite
        # then send warning that setup_link1d2d should be run again
        if overwrite_grid and "link1d2d" in self.data.data_vars:
            if grid_name == "mesh1d" or grid_name == "mesh2d":
                # TODO check if warning is enough or if we should remove to be sure?
                logger.warning(
                    f"{grid_name} grid was updated in self.mesh. "
                    "Re-run setup_link1d2d method to update the model 1D2D links."
                )

        # update related geoms if necessary: region - boundaries
        # the region is done in HydroMT Core
        if overwrite_grid or new_grid:
            # 1D boundaries
            if grid_name == "mesh1d":
                mesh1d_geom = workflows.get_boundaries_with_nodeid(
                    self.model.branches,
                    mesh_utils.network1d_nodes_geodataframe(
                        self.mesh_datasets["network1d"]
                    ),
                )
                self.model.geoms.set(mesh1d_geom, "boundaries")

    def set_link1d2d(
        self,
        link1d2d: xr.Dataset,
    ):
        """
        Add or replace the link1d2d in the model mesh.

        Parameters
        ----------
        link1d2d: xr.Dataset
            link1d2d dataset with variables: [link1d2d, link1d2d_ids,
            link1d2d_long_names, link1d2d_contact_type]
        """
        # Check if link1d2d already in self.mesh
        # FIXME current implementation of below does not support updating partial
        # 1d2d links. Either document or adapt. #1
        if "link1d2d" in self.data.data_vars:
            logger.info("Overwriting existing link1d2d in self.mesh.")
            self._data = self._data.drop_vars(
                [
                    "link1d2d",
                    "link1d2d_id",
                    "link1d2d_long_name",
                    "link1d2d_contact_type",
                ]
            )

        # Add link1d2d to mesh
        self._data = self._data.merge(link1d2d)
