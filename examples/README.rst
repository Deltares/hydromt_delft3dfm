This folder contains several iPython notebook examples for **HydroMT-Delft3D FM**. 

These examples can be run on your local machine. 

Local installation
------------------

To run these examples on your local machine you need a copy of the examples folder 
of the repository and an installation of HydroMT-Delft3D FM including some additional 
packages required to run the notebooks. 

1 - Install HydroMT-Delft3D FM 
******************************

The next step is to install all the python dependencies required to run the notebooks, including HydroMT and HydroMT-Delft3D FM.

**If you do not have hydromt-delft3dfm yet installed**, first create a new empty environment with the base hydromt installation:

.. code-block:: console

  $ conda create -n hydromt-delft3dfm -c conda-forge hydromt

Alternatively, you can also use mamba. Checkout more details in the `installation guide. <https://deltares.github.io/hydromt_delft3dfm/latest/getting_started/installation>`_

To run the notebooks, you need to install the ``examples`` version of HydroMT Delft3D FM using pip. The examples version installs additional dependencies
such as jupyter notebook to run the notebooks, matplotlib to plot etc. It is a more complete
installation of hydromt_delft3dfm. To install or update in an existing environment (example hydromt-delft3dfm environment), do:

.. code-block:: console

  $ conda activate hydromt-delft3dfm
  $ pip install "hydromt_delft3dfm[examples]"

2 - Download the content of the examples and notebooks
******************************************************
To run the examples locally, you will need to download the content of the hydromt_delft3dfm repository.
You have two options:

  1. Download and unzip the examples manually
  2. Clone the hydromt_delft3dfm GitHub repository

.. warning::

  Depending on your installed version of hydromt and hydromt_delft3dfm, you will need to download the correct versions of the examples.
  To check the version of hydromt_delft3dfm that you have installed, do:

  .. code-block:: console

    $ hydromt --models

    model plugins:
     - dflowfm (hydromt_delft3dfm 0.1.2)
    generic models (hydromt 0.9.0)
     - grid_model
     - vector_model
     - mesh_model
     - network_model

In the examples above, we see version 0.1.2 of hydromt_delft3dfm is installed and version 0.9.0 of hydromt.

**Option 1: manual download and unzip**

To manually download the examples on Windows, do (!replace with your own hydromt_delft3dfm version!):

.. code-block:: console

  $ curl https://github.com/Deltares/hydromt_delft3dfm/archive/refs/tags/v0.1.2.zip -O -L
  $ tar -xf v0.1.2.zip
  $ ren hydromt_delft3dfm-0.1.2 hydromt_delft3dfm

You can also download, unzip and rename manually if you prefer, rather than using the windows command prompt.

**Option 2: cloning the hydromt_delft3dfm repository**

For git users, you can also get the examples by cloning the hydromt_delft3dfm github repository and checking the version
you have installed:

.. code-block:: console

  $ git clone https://github.com/Deltares/hydromt_delft3dfm.git
  $ git checkout v0.1.2

3 - Running the examples
************************
Finally, start a jupyter lab server inside the **examples** folder 
after activating the **hydromt-delft3dfm** environment, see below.

Alternatively, you can run the notebooks from `Visual Studio Code <https://code.visualstudio.com/download>`_.

.. code-block:: console

  $ conda activate hydromt-delft3dfm
  $ cd hydromt_delft3dfm/examples
  $ jupyter lab
