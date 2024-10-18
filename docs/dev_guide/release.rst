.. _dev_release:

Releasing
=========

Release steps:

* Change version in __init__.py
* update the changelog
* Check dependencies in pyproj.toml
* run local testbank
* local check with: ``flit build`` and ``twine check dist/*``
* create PR and merge back to main
* create a new release: https://github.com/Deltares/hydromt_delft3dfm/releases/new)
* click ``choose a tag`` and type v+versionnumber (e.g. ``v0.1.0``), click ``create new tag: v0.1.0 on publish``
* set the release title to the tagname (e.g. ``v0.1.0``)
* click `Generate release notes` and replace the `What's Changed` info by a tagged link to ``docs/whats-new.md``
* if all is set, click ``Publish release``
* a release is created and the github action publishes it [on PyPI](https://pypi.org/project/dfm-tools)

Post release:

* Update the version to 0.1.1 (bump patch)
* Initialise changelog unrelease/added/fixed/changed
* Push and merge to main
