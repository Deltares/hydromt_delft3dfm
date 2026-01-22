.. currentmodule:: hydromt_delft3dfm

.. _api_reference:

#############
API reference
#############

.. _api_model:

DFlowFM model class
===================

Initialize
----------

.. autosummary::
   :toctree: _generated

   DFlowFM1D2DModel

.. _components:

Setup components
----------------

.. autosummary::
   :toctree: _generated

   DFlowFM1D2DModel.setup_config
   DFlowFM1D2DModel.setup_channels
   DFlowFM1D2DModel.setup_rivers_from_dem
   DFlowFM1D2DModel.setup_rivers
   DFlowFM1D2DModel.setup_pipes
   DFlowFM1D2DModel.setup_manholes
   DFlowFM1D2DModel.setup_1dboundary
   DFlowFM1D2DModel.setup_1dlateral_from_points
   DFlowFM1D2DModel.setup_1dlateral_from_polygons
   DFlowFM1D2DModel.setup_bridges
   DFlowFM1D2DModel.setup_culverts
   DFlowFM1D2DModel.setup_mesh2d
   DFlowFM1D2DModel.setup_mesh2d_refine
   DFlowFM1D2DModel.setup_link1d2d
   DFlowFM1D2DModel.setup_maps_from_rasterdataset
   DFlowFM1D2DModel.setup_maps_from_raster_reclass
   DFlowFM1D2DModel.setup_2dboundary
   DFlowFM1D2DModel.setup_rainfall_from_constant
   DFlowFM1D2DModel.setup_rainfall_from_uniform_timeseries

Attributes
----------

.. autosummary::
   :toctree: _generated

   DFlowFM1D2DModel.region
   DFlowFM1D2DModel.crs
   DFlowFM1D2DModel.bounds
   DFlowFM1D2DModel.res
   DFlowFM1D2DModel.root
   DFlowFM1D2DModel.config
   DFlowFM1D2DModel.maps
   DFlowFM1D2DModel.geoms
   DFlowFM1D2DModel.forcing
   DFlowFM1D2DModel.states
   DFlowFM1D2DModel.results
   DFlowFM1D2DModel.mesh
   DFlowFM1D2DModel.dfmmodel
   DFlowFM1D2DModel.dimr
   DFlowFM1D2DModel.branches
   DFlowFM1D2DModel.rivers
   DFlowFM1D2DModel.channels
   DFlowFM1D2DModel.pipes
   DFlowFM1D2DModel.opensystem
   DFlowFM1D2DModel.closedsystem
   DFlowFM1D2DModel.mesh_names
   DFlowFM1D2DModel.mesh_grids
   DFlowFM1D2DModel.mesh_datasets
   DFlowFM1D2DModel.mesh_gdf


High level methods
------------------

.. autosummary::
   :toctree: _generated

   DFlowFM1D2DModel.read
   DFlowFM1D2DModel.write
   DFlowFM1D2DModel.build
   DFlowFM1D2DModel.update
   DFlowFM1D2DModel.set_root
   DFlowFM1D2DModel.write_data_catalog

General methods
---------------

.. autosummary::
   :toctree: _generated

   DFlowFM1D2DModel.get_config
   DFlowFM1D2DModel.set_config
   DFlowFM1D2DModel.read_config
   DFlowFM1D2DModel.write_config

   DFlowFM1D2DModel.set_maps
   DFlowFM1D2DModel.read_maps
   DFlowFM1D2DModel.write_maps

   DFlowFM1D2DModel.set_geoms
   DFlowFM1D2DModel.read_geoms
   DFlowFM1D2DModel.write_geoms

   DFlowFM1D2DModel.set_forcing
   DFlowFM1D2DModel.read_forcing
   DFlowFM1D2DModel.write_forcing

   DFlowFM1D2DModel.set_states
   DFlowFM1D2DModel.read_states
   DFlowFM1D2DModel.write_states

   DFlowFM1D2DModel.set_results
   DFlowFM1D2DModel.read_results

   DFlowFM1D2DModel.get_mesh
   DFlowFM1D2DModel.set_mesh
   DFlowFM1D2DModel.set_link1d2d
   DFlowFM1D2DModel.read_mesh
   DFlowFM1D2DModel.write_mesh

   DFlowFM1D2DModel.set_branches

   DFlowFM1D2DModel.read_dimr
   DFlowFM1D2DModel.write_dimr

   DFlowFM1D2DModel.init_dfmmodel

   DFlowFM1D2DModel.get_model_time


.. _workflows:

DFlowFM1D2DModel workflows
======================

Boundaries
----------

.. autosummary::
   :toctree: _generated

   workflows.get_boundaries_with_nodeid
   workflows.select_boundary_type
   workflows.validate_boundaries
   workflows.compute_boundary_values
   workflows.compute_2dboundary_values
   workflows.compute_meteo_forcings
   workflows.compute_forcing_values_points
   workflows.compute_forcing_values_polygon
   workflows.get_geometry_coords_for_polygons

Branches
--------

.. autosummary::
   :toctree: _generated

   workflows.prepare_branches
   workflows.process_branches
   workflows.validate_branches
   workflows.add_branches
   workflows.find_nearest_branch
   workflows.update_data_columns_attributes
   workflows.update_data_columns_attribute_from_query
   workflows.snap_newbranches_to_branches_at_snappednodes
   workflows.snap_geom_to_branches_and_drop_nonsnapped

Crosssections
-------------

.. autosummary::
   :toctree: _generated

   workflows.prepare_default_friction_and_crosssection
   workflows.init_crosssections_options
   workflows.set_branch_crosssections
   workflows.set_xyz_crosssections
   workflows.set_point_crosssections
   workflows.add_crosssections

DEM
---

.. autosummary::
   :toctree: _generated

   workflows.invert_levels_from_dem
   workflows.get_river_bathymetry

Manholes
--------

.. autosummary::
   :toctree: _generated

   workflows.generate_manholes_on_branches

Mesh
----

.. autosummary::
   :toctree: _generated

   workflows.mesh1d_network1d_from_branches
   workflows.mesh1d_add_branch
   workflows.mesh2d_refine
   workflows.links1d2d_add_links_1d_to_2d
   workflows.links1d2d_add_links_2d_to_1d_embedded
   workflows.links1d2d_add_links_2d_to_1d_lateral

Region
------

.. autosummary::
   :toctree: _generated

   workflows.parse_region_geometry

Roughness
---------

.. autosummary::
   :toctree: _generated

   workflows.generate_roughness

Structures
----------

.. autosummary::
   :toctree: _generated

   workflows.prepare_1dstructures

.. _methods:

DFlowFM low-level methods
=========================

Input/Output methods
---------------------

.. autosummary::
   :toctree: _generated

   utils.read_branches_gui
   utils.write_branches_gui
   utils.read_crosssections
   utils.write_crosssections
   utils.read_friction
   utils.write_friction
   utils.read_structures
   utils.write_structures
   utils.read_manholes
   utils.write_manholes
   utils.read_1dboundary
   utils.write_1dboundary
   utils.read_1dlateral
   utils.write_1dlateral
   utils.read_2dboundary
   utils.write_2dboundary
   utils.read_meteo
   utils.write_meteo

Mesh conversion methods
-----------------------

.. autosummary::
   :toctree: _generated

   mesh_utils.hydrolib_network_from_mesh
   mesh_utils.mesh1d_network1d_from_hydrolib_network
   mesh_utils.links1d2d_from_hydrolib_network
   mesh_utils.mesh_from_hydrolib_network
   mesh_utils.mesh1d_nodes_geodataframe
   mesh_utils.network1d_nodes_geodataframe

Graph methods
-------------

.. autosummary::
   :toctree: _generated

   graph_utils.gpd_to_digraph
   graph_utils.get_endnodes_from_lines

GIS methods
-----------

.. autosummary::
   :toctree: _generated

   gis_utils.split_lines
   gis_utils.cut_pieces
   gis_utils.check_gpd_attributes
   gis_utils.update_data_columns_attributes_based_on_filter
   gis_utils.get_gdf_from_branches
