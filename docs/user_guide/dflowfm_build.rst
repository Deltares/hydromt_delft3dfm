.. _dflowfm_build:

Building a model
================

This plugin allows to build a complete model from available data. Once the configuration and
data libraries are set, you can build a model by using:

.. code-block:: console

    hydromt build dflowfm path/to/built_model -i dflowfm_build.yml -d data_sources.yml -vvv

.. _model_config:

Configuration file
------------------
Settings to build or update a Delft3DFM model are managed in a configuration file. In this file,
every option from each :ref:`model method <model_methods>` can be changed by the user
in its corresponding section.

Note that the order in which the components are listed in the yml file is important:


- When setting up a 1D model, one of the `setup_rivers`, `setup_channels` and `setup_pipes` should always be run first to determine the model topology.
- When setting up a 2D model, `setup_mesh2d` should always be run first.
- When setting up a 1D2D model, both of the above should be run first, before calling `setup_link1d2d`.


Below is an example yml file that can be used to build a complete Delft3DFM model
:download:`.yml file <../_examples/delft3dfm_build.yml>`. Each section corresponds
to a model component with the same name.

.. literalinclude:: ../_examples/delft3dfm_build.yml
   :language: yaml

Selecting data
--------------
Data sources in HydroMT are provided in one of several yaml libraries. These libraries contain required
information on the different data sources so that HydroMT can process them for the different models. There
are three ways for the user to select which data libraries to use:

- If no yaml file is selected, HydroMT will use the data stored in the
  `hydromt-artifacts <https://github.com/DirkEilander/hydromt-artifacts>`_
  which contains an extract of global data for a small region around the Piave river in Northern Italy.
- Another options for Deltares users is to select the deltares-data library (requires access to the Deltares
  P-drive). In the command lines examples below, this is done by adding either **--dd** or **--deltares-data**
  to the build / update command line.
- Finally, the user can prepare its own yaml libary (or libraries) (see
  `HydroMT documentation <https://deltares.github.io/hydromt/latest/_examples/prep_data_catalog.html>`_ to check the guidelines).
  These user libraries can be added either in the command line using the **-d** option and path/to/yaml or in the **yml/ini file**
  with the **data_libs** option in the [global] sections.

.. toctree::
    :hidden:

    Example: Build Delft3DFM 1D2D model <../_examples/build_1d2dmodel.ipynb>
