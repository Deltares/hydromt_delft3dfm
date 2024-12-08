==========
What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

UNRELEASED
==========

Added
-----
- Included new example notebook for basic HydroMT use through CLI. (PR #202)

Changed
-------

Fixed
-----

v0.3.0 (18 October 2024)
==========
This version includes many code updates to accomodate for newer versions of dependencies.
The only upper bound left is hydromt<1.

Added
-----
- Setup 1D laterals at branches from points and polygons. (PR #81)
- Added tests for Mesh workflows. (PR #133)
- Added tests for xyz cross-section workflows (PR #153).

Changed
-------
- Change default spacing in setup_channels from ``None`` to ``np.inf``. (PR #133)
- Added ``maxdist`` variable to setup_rivers and setup_channels. (PR #153)
- Renamed ``manhole_defaults_fn`` to ``manholes_defaults_fn`` in ``setup_manholes`` for consistency. (PR #187)
- No data values in maps from ``setup_maps_from_rasterdataset`` are now handled as -999.0. (PR #161)

Fixed
-----
- Bugfixing of reading of frictions (global), crosssections and boundaries when update. (PR #81)
- Support for xugrid>=0.9.0, meshkernel>=4.3.0, hydromt>=0.10,<1, pandas>=2. (PR #129)
- Fixing setup_links1d2d for 2d to 1d direction. (PR #133)
- Support for hydrolib-core>=0.8.0. (PR #139)
- Add support for Python 3.12. (PR #149)
- Fix writing of structures with newer (geo)pandas versions. (PR #151)
- Several bugfixes related to processing of cross-sections (PR #153)
- Support for geopandas v1 (PR #158)
- Support for latest version hydromt artifact data. (PR #160)
- Avoid sediment section in mdu so generated models can run in Delft3D FM Suite 2024.03 1D2D. (PR #184)
- fixed typo so ``setup_pipes()`` now allows field ``invlev_dn``. (PR #193)

v0.2.0 (20 November 2023)
=========================
Major dependency upgrade and add support for compound structures.

Added
-----
- Support reading and writing of compound structures. (PR#113)

Changed
-------
- Upgraded meshkernel dependency to version 3.0.0. (PR#109)
- Upgraded xugrid dependency to version 0.7.1. (PR#109)
- Upgraded hydrolib-core dependency to version 0.6.0. (PR#109)
- Support multiple structures at the same location. (PR#113)

v0.1.2 (20 October 2023)
========================
Added examples, updated documentation and Upgraded dependency.

Added
-----
- Basic examples of building and updateing 1D, 2D and 1D2D models. (PR#85)

Changed
-------
- Upgraded hydromt dependency to version 0.9.0. (PR#100)
- Updated documentation. (PR #97)

v0.1.1 (13 October 2023)
========================
Dependencies upgrade.

Changed
-------
- Upgraded meshkernel dependency to version 2.1.0. (PR#94)

v0.1.0 (22 September 2023)
==========================
First experimental release of functionalities to build and update Delft3D FM 1D2D models.

Added
-----
- Setup 1D network including rivers, channels and pipes;
- Setup 1D network components including manholes, bridges and culverts
- Setup 1D forcings including boundary and laterals.
- Setup 2D mesh and refine 2D mesh.
- Setup 2D maps including bedlevels, roughness, infiltration capacity, etc.
- Setup 2D forcings including boundary and precipitation.
- Setup 1D2D links.

.. _Keep a Changelog: http://keepachangelog.com/en/1.0.0/
.. _Semantic Versioning: http://semver.org/spec/v2.0.0.html
