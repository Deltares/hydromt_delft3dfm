.. _installation_guide:

==================
Installation Guide
==================

Prerequisites
=============
For more information about the prerequisites for an installation of the HydroMT package
and related dependencies, please visit the documentation of
`HydroMT core <https://deltares.github.io/hydromt/latest/getting_started/installation.html#installation-guide>`_

Compared to HydroMT, HydroMT Delft3D FM has additional dependencies, namely:

- `xugrid <https://github.com/Deltares/xugrid>`_
- `hydrolib-core <https://github.com/Deltares/HYDROLIB-core>`_
- `meshkernel-py <https://github.com/Deltares/MeshKernelPy>`_
- `networkx <https://networkx.org/>`_

Installation
============

HydroMT Delft3D FM is available on pypi.

Install HydroMT Delft3D FM in a new environment
----------------------------------------------

You can install HydroMT Delft3D FM in a new environment (recommended!) called `hydromt-delft3dfm`.
HydroMT Delft3D FM is available on pypi but not yet available on conda-forge.
Therefore, we recommend creating a conda environment and install everything with pip.

.. code-block:: console

  $ conda create -n hydromt-delft3dfm python=3.11 -c conda-forge

Then, activate the environment:

.. code-block:: console

  $ conda activate hydromt-delft3dfm

Finally, install HydroMT Delft3D FM from pypi using pip.

.. code-block:: console

  $ pip install hydromt-delft3dfm

.. Tip::

    If you already have this environment with this name either remove it with
    `conda remove -n hydromt-delft3dfm --all` **or** set a new name for the environment
    by adding `-n <name>` to the line below.


Developer install
==================
To be able to test and develop the HydroMT-delft3dfm package see instructions in the :ref:`Developer installation guide <dev_env>`.
