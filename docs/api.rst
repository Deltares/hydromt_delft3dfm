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
   DFlowFMModel.bounds
   DFlowFMModel.crs
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
   DFlowFMModel.crosssections
   DFlowFMModel.boundaries

High level methods
------------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.read
   DFlowFMModel.write
   DFlowFMModel.build
   DFlowFMModel.update
   DFlowFMModel.set_root
   DFlowFMModel._model_has_2d
   DFlowFMModel._model_has_1d

General methods
---------------

.. autosummary::
   :toctree: _generated


   DFlowFMModel.setup_config
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

   DFlowFMModel.set_mesh
   DFlowFMModel.set_link1d2d
   DFlowFMModel.read_mesh
   DFlowFMModel.write_mesh

   DFlowFMModel.read_dimr
   DFlowFMModel.write_dimr

   DFlowFMModel.set_branches
   DFlowFMModel.add_branches

   DFlowFMModel.add_crosssections

   DFlowFMModel.get_boundaries
   DFlowFMModel.set_boundaries

   DFlowFMModel.get_model_time


.. _workflows:

DFlowFMModel workflows
======================

.. autosummary::
   :toctree: _generated

   workflows.branches


.. _methods:

DFlowFM low-level methods
=========================

Input/Output methods
---------------------

.. autosummary::
   :toctree: _generated

   utils.read_branches_gui

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
