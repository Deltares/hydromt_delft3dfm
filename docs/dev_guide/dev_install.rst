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

Then, make and activate a new hydromt-delft3dfm conda environment based on the envs/hydromt-delft3dfm.yml
file contained in the repository:

.. code-block:: console

    $ conda env create -f envs/hydromt-delft3dfm.yml
    $ conda activate hydromt-delft3dfm

If you wish to make changes in HydroMT-Delft3DFM, you should make an editable install of the plugin.
This can be done with:

.. code-block:: console

    $ pip install -e .

If you encounter issues with the installation of some packages, you might consider cleaning conda to remove unused packages and caches.
This can be done through the following command from your base environment:

.. code-block:: console

    $ conda clean -a
