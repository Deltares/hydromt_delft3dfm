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

   DFlowFMModel

.. _components:

Setup components
----------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.setup_config
   DFlowFMModel.setup_channels
   DFlowFMModel.setup_rivers_from_dem
   DFlowFMModel.setup_rivers
   DFlowFMModel.setup_pipes
   DFlowFMModel.setup_manholes
   DFlowFMModel.setup_1dboundary
   DFlowFMModel.setup_1dlateral_from_points
   DFlowFMModel.setup_1dlateral_from_polygons
   DFlowFMModel.setup_bridges
   DFlowFMModel.setup_culverts
   DFlowFMModel.setup_mesh2d
   DFlowFMModel.setup_mesh2d_refine
   DFlowFMModel.setup_link1d2d
   DFlowFMModel.setup_maps_from_rasterdataset
   DFlowFMModel.setup_maps_from_raster_reclass
   DFlowFMModel.setup_2dboundary
   DFlowFMModel.setup_rainfall_from_constant
   DFlowFMModel.setup_rainfall_from_uniform_timeseries

Attributes
----------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.region
   DFlowFMModel.crs
   DFlowFMModel.bounds
   DFlowFMModel.root
   DFlowFMModel.mdu
   DFlowFMModel.inifield
   DFlowFMModel.geoms
   DFlowFMModel.forcing
   DFlowFMModel.mesh
   DFlowFMModel.dfmmodel
   DFlowFMModel.dimr
   DFlowFMModel.branches
   DFlowFMModel.rivers
   DFlowFMModel.channels
   DFlowFMModel.pipes
   DFlowFMModel.opensystem
   DFlowFMModel.closedsystem

High level methods
------------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.read
   DFlowFMModel.write
   DFlowFMModel.build
   DFlowFMModel.update
   DFlowFMModel.write_data_catalog

General methods
---------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.set_branches

   DFlowFMModel.init_dfmmodel

   DFlowFMModel.get_model_time


.. _workflows:

DFlowFMModel workflows
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

   utils.io_utils.read_branches_gui
   utils.io_utils.write_branches_gui
   utils.io_utils.read_crosssections
   utils.io_utils.write_crosssections
   utils.io_utils.read_friction
   utils.io_utils.write_friction
   utils.io_utils.read_structures
   utils.io_utils.write_structures
   utils.io_utils.read_manholes
   utils.io_utils.write_manholes
   utils.io_utils.read_1dboundary
   utils.io_utils.write_1dboundary
   utils.io_utils.read_1dlateral
   utils.io_utils.write_1dlateral
   utils.io_utils.read_2dboundary
   utils.io_utils.write_2dboundary
   utils.io_utils.read_meteo
   utils.io_utils.write_meteo

Mesh conversion methods
-----------------------

.. autosummary::
   :toctree: _generated

   utils.mesh_utils.hydrolib_network_from_mesh
   utils.mesh_utils.mesh1d_network1d_from_hydrolib_network
   utils.mesh_utils.links1d2d_from_hydrolib_network
   utils.mesh_utils.mesh_from_hydrolib_network
   utils.mesh_utils.mesh1d_nodes_geodataframe
   utils.mesh_utils.network1d_nodes_geodataframe

Graph methods
-------------

.. autosummary::
   :toctree: _generated

   utils.graph_utils.gpd_to_digraph
   utils.graph_utils.get_endnodes_from_lines

GIS methods
-----------

.. autosummary::
   :toctree: _generated

   utils.gis_utils.split_lines
   utils.gis_utils.cut_pieces
   utils.gis_utils.check_gpd_attributes
   utils.gis_utils.update_data_columns_attributes_based_on_filter
   utils.gis_utils.get_gdf_from_branches
