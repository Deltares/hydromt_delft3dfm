.. _readme:

=================================================
HydroMT-Delft3D FM: Delft3D FM plugin for HydroMT
=================================================

|ci| |codecov| |black| |docs_latest|  |license| |sonarcloud|

What is the HydroMT-Delft3D FM plugin
-------------------------------------
HydroMT_ (Hydro Model Tools) is an open-source Python package that facilitates the process of
building and analyzing spatial geoscientific models with a focus on water system models.
It does so by automating the workflow to go from raw data to a complete model instance which
is ready to run and to analyze model results once the simulation has finished.
This plugin provides an implementation of the model API for the Delft3D FM_ model.

Why HydroMT-Delft3D FM?
-----------------------
Setting up distributed hydrological models typically requires many (manual) steps
to process input data and might therefore be time consuming and hard to reproduce.
Especially improving models based on global-local geospatial datasets, which are
rapidly becoming available at increasingly high resolutions, might be challenging.
HydroMT-Delft3D FM aims to make the Delft3D FM model building and updating processes **fast**, **modular** and **reproducible**
and to facilitate the analysis of the model results.

How to use HydroMT-Delft3D FM?
------------------------------
The HydroMT-Delft3D FM plugin can be used as a **command line** application, which provides commands to *build*,
and *update* a Delft3D FM model with a single line, or **from python** to exploit its rich interface.
You can learn more about how to use HydroMT-Delft3D FM in its `online documentation. <https://deltares.github.io/hydromt_delft3dfm/latest/getting_started/intro>`_
For a smooth installing experience we recommend installing HydroMT-Delft3D FM and its dependencies
from conda-forge in a clean environment, see `installation guide. <https://deltares.github.io/hydromt_delft3dfm/latest/getting_started/installation>`_

How to cite?
------------
For publications, please cite our work using the DOI provided in the Zenodo badge that points to the latest release.

How to contribute?
-------------------
If you find any issues in the code or documentation feel free to leave an issue on the `github issue tracker. <https://github.com/Deltares/hydromt_delft3dfm/issues>`_
You can find information about how to contribute to the HydroMT project at our `contributing page. <https://deltares.github.io/hydromt/latest/dev/contributing>`_

HydroMT seeks active contribution from the (hydro) geoscientific community.
So far, it has been developed and tested with a range of `Deltares <https://www.deltares.nl/en/>`_ models, but
we believe it is applicable to a much wider set of geoscientific models and are
happy to discuss how it can be implemented for your model.

.. _Hydromt: https://deltares.github.io/hydromt/latest/
.. _FM: https://oss.deltares.nl/web/delft3dfm

.. |ci| image:: https://github.com/Deltares/hydromt_delft3dfm/actions/workflows/ci.yml/badge.svg?branch=main
    :alt: ci
    :target: https://github.com/Deltares/hydromt_delft3dfm/actions/workflows/ci.yml

.. |codecov| image:: https://codecov.io/gh/Deltares/hydromt_delft3dfm/branch/main/graph/badge.svg?token=ss3EgmwHhH
    :target: https://codecov.io/gh/Deltares/hydromt_delft3dfm

.. |black|  image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :alt: Formatter
    :target: https://github.com/psf/black

.. |docs_latest| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
    :target: https://deltares.github.io/hydromt_delft3dfm/latest
    :alt: Latest developers docs

.. |license| image:: https://img.shields.io/github/license/Deltares/hydromt_delft3dfm
    :alt: License
    :target: https://github.com/Deltares/hydromt_delft3dfm/blob/main/LICENSE

.. |sonarcloud| image:: https://sonarcloud.io/api/project_badges/measure?project=Deltares_hydromt_delft3dfm&metric=alert_status
    :alt: Quality Gate Status
    :target: https://sonarcloud.io/summary/new_code?id=Deltares_hydromt_delft3dfm
