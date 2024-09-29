==========
What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

Unreleased
==========
This version includes many code updates to accomodate for newer versions of dependencies.
This makes the code more robust and future-proof.

Added
-----
- Setup 1D laterals at branches from points and polygons. (PR #81)
- Adding tests for Mesh workflows. (PR #133)

Changed
-------
- Change default spacing in setup_channels from ``None`` to ``np.inf``. (PR #133)
- No data values in maps from ``setup_maps_from_rasterdataset`` are now handled as -999.0. (PR #161)

Fixed
-----
- Bugfixing of reading of frictions (global), crosssections and boundaries when update. (PR #81)
- Fixing bug related to changes to pandas TimeDelta formatting, see also https://pandas.pydata.org/docs/whatsnew/v2.2.0.html#other-deprecations. (PR #129)
- Fixing setup_links1d2d for 2d to 1d direction. (PR #133)

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
