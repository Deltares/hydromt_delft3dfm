"""Custom Delft3D FM geoms component."""

# TODO: channels, manholes etc should become their own components
import itertools
import logging
from os.path import dirname, join

import numpy as np
import pandas as pd
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import GeomsComponent, ModelComponent

from hydromt_delft3dfm.utils import io_utils, mesh_utils

__all__ = ["Delft3DFMGeomsComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class Delft3DFMGeomsComponent(GeomsComponent):
    """
    Manage the Delft3D FM geoms files for model geometries.

    It extends the base GeomsComponent from hydromt and consists of a dictionary of
    geopandas GeoDataFrames.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "geoms/{name}.geojson",
        region_component: str | None = None,
    ):
        """Initialize Delft3DFMGeomsComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            by default "geoms/{name}.geojson", i.e. one file per geodataframe in
            the data dictionary.
        region_component : str, optional
            The name of the region component to use as reference for this component's
            region. If None, the region will be set to the union of all geometries in
            the data dictionary.
        """
        super().__init__(
            model=model,
            filename=filename,
            region_component=region_component,
        )

    ### I/O methods ###
    @hydromt_step
    def read(self) -> None:  # FIXME: gives an error when only 2D model.
        """
        Read model geometries files at <root>/<geoms> and add to geoms property.

        For branches / boundaries etc... the reading of hydrolib-core objects happens
        in read_mesh. There the geoms geojson copies are re-set based on dflowfm files
        content.
        """
        self.root.is_reading_mode()
        super().read()

        if self.model.dfmmodel.geometry.crosslocfile is not None:
            # Read cross-sections and friction
            # Add crosssections properties, should be done before friction
            # Branches are needed do derive locations,
            # self.branches should start the read if not done yet
            logger.info("Reading cross-sections files")
            crosssections = io_utils.read_crosssections(
                self.model.branches, self.model.dfmmodel
            )

            # Add friction properties from roughness files
            # logger.info("Reading friction files")
            crosssections = io_utils.read_friction(crosssections, self.model.dfmmodel)
            self.set(crosssections, "crosssections")

        # Read manholes
        if self.model.dfmmodel.geometry.storagenodefile is not None:
            logger.info("Reading manholes file")
            network1d_nodes = mesh_utils.network1d_nodes_geodataframe(
                self.model.mesh.mesh_datasets["network1d"]
            )
            manholes = io_utils.read_manholes(network1d_nodes, self.model.dfmmodel)
            self.set(manholes, "manholes")

        # Read structures
        if self.model.dfmmodel.geometry.structurefile is not None:
            logger.info("Reading structures file")
            structures = io_utils.read_structures(
                self.model.branches, self.model.dfmmodel
            )
            for st_type in structures["type"].unique():
                self.set(structures[structures["type"] == st_type], f"{st_type}s")

    @hydromt_step
    def write(self, write_mesh_gdf: bool = True) -> None:
        """Write model geometries to a GeoJSON file at <root>/<geoms>."""
        self.root.is_writing_mode()

        # Optional: also write mesh_gdf object
        if write_mesh_gdf:
            for name, gdf in self.model.mesh.mesh_gdf.items():
                self.set(gdf, name)

        if len(self.data) == 0:
            logger.info(
                f"{self.model.name}.{self.name_in_model}: "
                "No geoms data found, skip writing."
            )
            return

        # Write geojson equivalent of all objects.
        # Note that these files are not directly used when updating the model
        super().write(precision=None)

        # Write dfm files
        savedir = dirname(join(self.root.path, self.model.mdu._filename))

        # Write cross-sections (inc. friction)
        if "crosssections" in self._data:
            # Crosssections
            gdf_crs = self.data["crosssections"]
            logger.info("Writting cross-sections files crsdef and crsloc")
            crsdef_fn, crsloc_fn = io_utils.write_crosssections(gdf_crs, savedir)
            self.model.mdu.set("geometry.crossdeffile", crsdef_fn)
            self.model.mdu.set("geometry.crosslocfile", crsloc_fn)
            # Friction
            logger.info("Writting friction file(s)")
            friction_fns = io_utils.write_friction(gdf_crs, savedir)
            self.model.mdu.set("geometry.frictfile", ";".join(friction_fns))

        # Write structures
        # Manholes
        if "manholes" in self._data:
            logger.info("Writting manholes file.")
            storage_fn = io_utils.write_manholes(
                self.data["manholes"],
                savedir,
            )
            self.model.mdu.set("geometry.storagenodefile", storage_fn)

        # Write structures
        existing_structures = [st for st in ["bridges", "culverts"] if st in self.data]
        if len(existing_structures) > 0:
            # combine all structures
            structures = []
            for st in existing_structures:
                structures.append(self.data.get(st).to_dict("records"))
            structures = list(itertools.chain.from_iterable(structures))
            structures = pd.DataFrame(structures).replace(np.nan, None)
            # write
            logger.info("Writting structures file.")
            structures_fn = io_utils.write_structures(
                structures,
                savedir,
            )
            self.model.mdu.set("geometry.structurefile", structures_fn)

        # write hydromt
        # Optional: also write mesh_gdf object
        if write_mesh_gdf:
            for name, gdf in self.model.mesh.mesh_gdf.items():
                self.set(gdf, name)

        # Write geojson equivalent of all objects.
        # NOTE these files are not used for model update.
        # convert any list in geoms to strings
        def convert_lists_to_strings(df):
            for column_name in df.columns:
                if df[column_name].apply(isinstance, args=(list,)).any():
                    df[column_name] = df[column_name].apply(
                        lambda x: " ".join(f"{x}") if isinstance(x, list) else x
                    )
            return df

        for name in self.data:
            self.set(convert_lists_to_strings(self.data[name]), name)

    def test_equal(self, other: ModelComponent) -> tuple[bool, dict[str, str]]:
        """Test if two GeomsComponents are equal.

        Parameters
        ----------
        other: GeomsComponent
            The other GeomsComponent to compare with.

        Returns
        -------
        tuple[bool, dict[str, str]]
            True if the components are equal, and a dict with the associated errors per
            property checked.
        """
        eq, errors = super().test_equal(other)

        return eq, {"geoms": errors}
