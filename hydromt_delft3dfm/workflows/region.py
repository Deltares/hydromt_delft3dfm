"""Workflows to parse region for Delft3D FM model."""


from hydromt.model.processes.region import parse_region_bbox, parse_region_geom
from pyproj.crs import CRS

__all__ = [
    "parse_region_geometry",
]


def parse_region_geometry(
    region: dict,
    crs: CRS,
):
    """Parse hydromt region argument into region geometry."""
    kind = next(iter(region))
    if kind == "bbox":
        geom = parse_region_bbox(region)
    elif kind == "geom":
        geom = parse_region_geom(region)
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
