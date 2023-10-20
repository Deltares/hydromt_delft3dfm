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

- `xugrid <https://github.com/Deltares/xugrid>`_
- `hydrolib-core <https://github.com/Deltares/HYDROLIB-core>`_
- `meshkernel-py <https://github.com/Deltares/MeshKernelPy>`_
- `networkx <https://networkx.org/>`_

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
-------------------------------------------------------------

You can install HydroMT-Delft3DFM in a new environment called `hydromt-delft3dfm`.
HydroMT-Delft3DFM is not yet available on conda-forge so we recommend installling HydroMT (core) first
via conda-forge and then hydromt-delft3dfm and the additionnal libraries via pip.

.. code-block:: console

  $ conda env create -n hydromt-delft3dfm -c conda-forge hydromt

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


Developer install
==================
To be able to test and develop the HydroMT-delft3dfm package see instructions in the :ref:`Developer installation guide <dev_env>`.
