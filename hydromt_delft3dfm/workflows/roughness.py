"""Workflows to prepare roughness for Delft3D-FM model."""

import logging

import geopandas as gpd

logger = logging.getLogger("hydromt")


__all__ = ["generate_roughness"]


def generate_roughness(roughness: gpd.GeoDataFrame):
    """Generate roughness ID column based on frictiontype and frictionvalue."""
    roughness["frictionid"] = roughness.apply(
        lambda x: "%s_%s" % (x["frictiontype"], x["frictionvalue"]), axis=1
    )
    roughness_only = roughness[
        [roughness.index.name, "geometry", "frictiontype", "frictionvalue"]
    ]

    return roughness_only, roughness
