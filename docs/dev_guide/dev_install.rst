.. _dev_env:

Developer's environment
=======================
If you want to download the HydroMT-Delft3DFM plugin directly from git to easily have access to the latest developments or
make changes to the code you can use the following steps.

First, clone the HydroMT-Delft3DFM plugin ``git`` repo from
`github <https://github.com/Deltares/hydromt_delft3dfm>`_, then navigate into the
the code folder (where the envs folder and pyproject.toml are located):

.. code-block:: console

    $ git clone https://github.com/Deltares/hydromt_delft3dfm.git
    $ cd hydromt_delft3dfm

Then, create and activate a new conda environment and install hydromt_delft3dfm as developer:

.. code-block:: console

    $ conda create -n hydromt-delft3dfm python=3.11 -c conda-forge
    $ conda activate hydromt-delft3dfm
    $ pip install -e .
