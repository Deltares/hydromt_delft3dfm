.. _dev_env:

Developer's guide
=================
If you want to download the HydroMT-Delft3DFM plugin directly from git to easily have access to the latest developments or
make changes to the code you can use the following steps.

Install
-------

First, install Pixi: `<https://pixi.prefix.dev/latest/installation>`_

Then, clone the HydroMT-Delft3DFM plugin ``git`` repo from
`github <https://github.com/Deltares/hydromt_delft3dfm>`_, then navigate into the
the code folder (where the pyproject.toml is located):

.. code-block:: console

    $ git clone https://github.com/Deltares/hydromt_delft3dfm.git
    $ cd hydromt_delft3dfm

If git is not recognized as a command, first install Git from `<https://git-scm.com/install>`_. VSCode and PyCharm should be bundled with a git extension so a manual installation is not always necessary.

Then, create and activate a new pixi environment. This  includes a developers installation of hydromt_delft3dfm and several environments as specified in pyproject.toml:

.. code-block:: console

    $ pixi install

Test
----

Running the tests can be done within any of the environment specified in the pyproject.toml (defaults to `default`):

.. code-block:: console

    $ pixi run -e default test
    $ pixi run -e default pytest

Updating the lockfile
---------------------

If you add any dependencies or change anything in the package configuration, you have to update the lockfile. This is also done (if needed) when running other pixi commands:

.. code-block:: console

    $ pixi lock

Generating the docs
-------------------

Generating the docs is added as a pixi task, which executes sphinx-build:

.. code-block:: console

    $ pixi run docs-build

If any developer's information is missing, it might be documented in the `hydromt developer's guide <https://deltares.github.io/hydromt/stable/dev/core_dev/index.html>`_.
