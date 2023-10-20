==========
What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

Unreleased
==========

Added
-----
- Basic examples of building and updateing 1D, 2D and 1D2D models. (PR#85)

Changed
-------
- Upgraded hydromt dependency to version 0.9.0. (PR#100) 

Fixed
-----

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
