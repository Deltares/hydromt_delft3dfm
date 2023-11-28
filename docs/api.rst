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
   DFlowFMModel.res
   DFlowFMModel.root
   DFlowFMModel.config
   DFlowFMModel.maps
   DFlowFMModel.geoms
   DFlowFMModel.forcing
   DFlowFMModel.states
   DFlowFMModel.results
   DFlowFMModel.mesh
   DFlowFMModel.dfmmodel
   DFlowFMModel.dimr
   DFlowFMModel.branches
   DFlowFMModel.rivers
   DFlowFMModel.channels
   DFlowFMModel.pipes
   DFlowFMModel.opensystem
   DFlowFMModel.closedsystem
   DFlowFMModel.mesh_names
   DFlowFMModel.mesh_grids
   DFlowFMModel.mesh_datasets
   DFlowFMModel.mesh_gdf


High level methods
------------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.read
   DFlowFMModel.write
   DFlowFMModel.build
   DFlowFMModel.update
   DFlowFMModel.set_root
   DFlowFMModel.write_data_catalog

General methods
---------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.get_config
   DFlowFMModel.set_config
   DFlowFMModel.read_config
   DFlowFMModel.write_config

   DFlowFMModel.set_maps
   DFlowFMModel.read_maps
   DFlowFMModel.write_maps

   DFlowFMModel.set_geoms
   DFlowFMModel.read_geoms
   DFlowFMModel.write_geoms

   DFlowFMModel.set_forcing
   DFlowFMModel.read_forcing
   DFlowFMModel.write_forcing

   DFlowFMModel.set_states
   DFlowFMModel.read_states
   DFlowFMModel.write_states

   DFlowFMModel.set_results
   DFlowFMModel.read_results

   DFlowFMModel.get_mesh
   DFlowFMModel.set_mesh
   DFlowFMModel.set_link1d2d
   DFlowFMModel.read_mesh
   DFlowFMModel.write_mesh

   DFlowFMModel.set_branches

   DFlowFMModel.read_dimr
   DFlowFMModel.write_dimr

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
   workflows.generate_boundaries_from_branches
   workflows.select_boundary_type
   workflows.validate_boundaries
   workflows.compute_boundary_values
   workflows.compute_2dboundary_values
   workflows.compute_meteo_forcings

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
