.. currentmodule:: hydromt_delft3dfm

.. _model_set_up:

============================
Model methods and components
============================

The HydroMT-Delft3D FM plugin helps you preparing or updating several methods of a Delft3D FM model such as model topology, map data and forcing.
The main interactions are available from the HydroMT Command Line Interface and allow you to configure
HydroMT in order to build or update Delft3D FM models.

When building or updating a model from command line a model setup
:ref:`configuration <model_config>` (.yml file) with model methods and options and, optionally,
a `data sources <https://deltares.github.io/hydromt/latest/user_guide/data_prepare_cat.html>`_ (.yml) file should be prepared.

.. _model_methods:

Model setup methods
===================

An overview of the available Delft3D FM model setup methods
is provided in the table below. When using HydroMT from the command line only the
setup methods are exposed. Click on
a specific method see its documentation.

.. list-table::
    :widths: 20 55
    :header-rows: 1
    :stub-columns: 1

    * - Method
      - Explanation
    * - :py:func:`~DFlowFMModel.setup_config`
      - Update config with a dictionary
    * - :py:func:`~DFlowFMModel.setup_rivers`
      - This component sets the 1D river branches with parameters.
    * - :py:func:`~DFlowFMModel.setup_channels`
      - This component sets the 1D channel branches with parameters.
    * - :py:func:`~DFlowFMModel.setup_rivers_from_dem`
      - This component sets the 1D river branches with parameters derived from hydrography.
    * - :py:func:`~DFlowFMModel.setup_pipes`
      - This component sets the 1D pipe branches with parameters.
    * - :py:func:`~DFlowFMModel.setup_manholes`
      - This component adds manholes with parameters to 1D pipes.
    * - :py:func:`~DFlowFMModel.setup_bridges`
      - This component adds bridges with parameters to 1D branches.
    * - :py:func:`~DFlowFMModel.setup_culverts`
      - This component adds culverts with parameters to 1D branches.
    * - :py:func:`~DFlowFMModel.setup_mesh2d`
      - This component sets a 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_mesh2d_refine`
      - This component refines the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_link1d2d`
      - This component sets 1d2d links that link the 1D branchs to the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_maps_from_rasterdataset`
      - This component adds parameter maps to the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_maps_from_raster_reclass`
      - This component adds parameter maps that are derived by reclass existing maps to the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_1dboundary`
      -  Setup a 1D boundary forcing to the 1D branches.
    * - :py:func:`~DFlowFMModel.setup_2dboundary`
      -  Setup a 2D boundary forcing to the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_2dboundary_from_lines`
      -  Setup a 2D boundary forcing to the 2D mesh from lines.
    * - :py:func:`~DFlowFMModel.setup_rainfall_from_constant`
      -  Setup a constant precipitation forcing to the 2D mesh.
    * - :py:func:`~DFlowFMModel.setup_rainfall_from_uniform_timeseries`
      -  Setup a spatial uniform precipitation forcing to the 2D mesh.


.. _model_components:

Model components
================

The following table provides an overview of which :py:class:`~hydromt_delft3dfm.DFlowFMModel`
component contains which Delft3DFM in- and output files. The files are read and written with the associated
read- and write- methods, i.e. :py:func:`~DFlowFMModel.read_config`
and :py:func:`~DFlowFMModel.write_config` for the
:py:attr:`~DFlowFMModel.config` component.


.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - :py:class:`~hydromt_delft3dfm.DFlowFMModel` component
     - Delft3DFM files, saved in dflowfm folder.
   * - :py:attr:`~hydromt_delft3dfm.DFlowFMModel.config`
     - DFlowFM.mdu, saved in dflowfm folder.
   * - :py:attr:`~hydromt_delft3dfm.DFlowFMModel.geoms`
     - 1D model geometries, saved in the geoms folder.
   * - :py:attr:`~hydromt_delft3dfm.DFlowFMModel.maps`
     - 2D model data on regular grid, saved in the maps folder
   * - :py:attr:`~hydromt_delft3dfm.DFlowFMModel.mesh`
     - unstructured mesh schematisation and data, saved in the dflowfm folder
   * - :py:attr:`~hydromt_delft3dfm.DFlowFMModel.forcing`
     - 1D and 2D forcing, saved in dflowfm  folder.
