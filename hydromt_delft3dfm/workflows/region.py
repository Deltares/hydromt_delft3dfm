"""Workflows to parse region for Delft3D-FM model."""


import geopandas as gpd
from hydromt.model.processes.region import parse_region_geom
from pyproj.crs import CRS
from shapely.geometry import box

__all__ = [
    "parse_region_geometry",
]


def parse_region_geometry(
    region: dict,
    crs: CRS,
):
    """Parse hydromt stype region argument into region geometry."""
    kind, region = parse_region_geom(region)
    if kind == "bbox":
        bbox = region["bbox"]
        geom = gpd.GeoDataFrame(geometry=[box(*bbox)], crs=4326)
    elif kind == "geom":
        geom = region["geom"]
        if geom.crs is None:
            raise ValueError('Model region "geom" has no CRS')
    else:
        raise ValueError(
            f"Region for mesh1d must of kind [bbox, geom], kind {kind} "
            "not understood."
        )

    if geom.crs != crs:
        geom = geom.to_crs(crs)

    return geom
