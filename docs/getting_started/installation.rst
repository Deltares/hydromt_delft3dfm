.. _installation_guide:

==================
Installation Guide
==================

Prerequisites
=============
For more information about the prerequisites for an installation of the HydroMT package
and related dependencies, please visit the documentation of
`HydroMT core <https://deltares.github.io/hydromt/latest/getting_started/installation.html#installation-guide>`_

Compared to HydroMT, HydroMT-delft3dfm has additional dependencies, namely:

- `toml <https://github.com/uiri/toml>`_
- `pcraster <https://pcraster.geo.uu.nl>`_ (optional)
- `gwwapi <https://github.com/global-water-watch/gww-api>`_ (optional)
- `hydroengine <https://github.com/openearth/hydro-engine>`_ (optional)

If you already have a python & conda installation but do not yet have mamba installed,
we recommend installing it into your *base* environment using:

.. code-block:: console

  $ conda install mamba -n base -c conda-forge


Installation
============

HydroMT-Delft3DFM is available from pypi.
We recommend installing using mamba/conda from conda-forge in a new environment.
If conda is prefered, we recommond install libmamba as the solver. See how to `here <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_. 

.. Note::

    In the commands below you can exchange `mamba` for `conda`, see
    `here <https://deltares.github.io/hydromt/latest/getting_started/installation.html#installation-guide>`_
    for the difference between both.

Install HydroMT-Delft3DFM in a new environment (recommended!)
---------------------------------------------------------

You can install HydroMT-Delft3DFM in a new environment called `hydromt-delft3dfm` together with
all optional (see above) and a few additional dependencies with:

.. code-block:: console

  $ conda env create -f https://raw.githubusercontent.com/Deltares/hydromt_delft3dfm/54-examples/envs/hydromt-delft3dfm-min.yml

Then, activate the environment (as stated by mamba/conda) to start making use of HydroMT-delft3dfm:

.. code-block:: console

  conda activate hydromt-delft3dfm

Finally, install Hydromt-Delft3DFM using pypi.

.. code-block:: console

  pip install hydromt-delft3dfm

.. Tip::

    If you already have this environment with this name either remove it with
    `conda env remove -n hydromt-delft3dfm` **or** set a new name for the environment
    by adding `-n <name>` to the line below.

.. Install HydroMT-delft3dfm in an existing environment
.. ------------------------------------------------

.. To install HydroMT-Delft3DFM in an existing environment execute the command below
.. where you replace `<environment_name>` with the name of the existing environment.
.. Note that if some dependencies are not installed from conda-forge but from other
.. channels the installation may fail.

.. .. code-block:: console

..    $ conda activate <environment_name>
..    $ pip install hydromt_delft3dfm

.. .. Note::

..     Please take into account that hydromt is now installed from unreleased version on github.

.. .. code-block:: console

..   $ pip install git+https://github.com/Deltares/hydromt

Developer install
==================
To be able to test and develop the HydroMT-delft3dfm package see instructions in the :ref:`Developer installation guide <dev_env>`.
