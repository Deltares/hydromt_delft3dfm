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

Attributes
----------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.region
   DFlowFMModel.crs
   DFlowFMModel.res
   DFlowFMModel.root
   DFlowFMModel.config
   DFlowFMModel.staticmaps
   DFlowFMModel.staticgeoms
   DFlowFMModel.forcing
   DFlowFMModel.states
   DFlowFMModel.results


High level methods
------------------

.. autosummary::
   :toctree: _generated

   DFlowFMModel.read
   DFlowFMModel.write
   DFlowFMModel.build
   DFlowFMModel.update
   DFlowFMModel.set_root

General methods
---------------

.. autosummary::
   :toctree: _generated


   DFlowFMModel.setup_config
   DFlowFMModel.get_config
   DFlowFMModel.set_config
   DFlowFMModel.read_config
   DFlowFMModel.write_config

   DFlowFMModel.set_staticmaps
   DFlowFMModel.read_staticmaps
   DFlowFMModel.write_staticmaps
   DFlowFMModel.clip_staticmaps
   DFlowFMModel.set_staticgeoms
   DFlowFMModel.read_staticgeoms
   DFlowFMModel.write_staticgeoms

   DFlowFMModel.set_forcing
   DFlowFMModel.read_forcing
   DFlowFMModel.write_forcing
   DFlowFMModel.clip_forcing

   DFlowFMModel.set_states
   DFlowFMModel.read_states
   DFlowFMModel.write_states

   DFlowFMModel.set_results
   DFlowFMModel.read_results


.. _workflows:

Wflow workflows
===============

.. autosummary::
   :toctree: _generated

   workflows.branches


.. _methods:

DFlowFM low-level methods
=======================

Input/Output methods
---------------------

.. autosummary::
   :toctree: _generated

   utils.read_branches_gui
