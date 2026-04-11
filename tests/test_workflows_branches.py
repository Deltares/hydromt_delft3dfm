import pytest
import geopandas as gpd
from hydromt_delft3dfm.workflows.branches import reduce_gdf_precision
from shapely import LineString, Point


def test_reduce_gdf_precision():
    # test on a single branch with 3 cross-sections (id = 1 .. 3)
    linestrings = gpd.GeoDataFrame(
        data={"frictionid": ["Manning_0.023"], "frictiontype": ["Manning"], "frictionvalue": [0.023]}, 
        geometry=[LineString([[0, 0], [1000, 1000]])],
        crs=28992)
    linestrings.index = linestrings.index.astype(str)

    points = gpd.GeoDataFrame(
        geometry=[Point(-10, 1), Point(0, 1), Point(10, 2), Point(20, 2)],
        data={"crsid": [1, 1, 1, 1], "order": [1, 2, 3, 4], "z": [10, 0, 1, 11]}, crs=28992)
    points.index = points.index.astype(str)

    mixed = gpd.GeoDataFrame(
        geometry=[Point(-10, 1), Point(0, 1), Point(10, 2), LineString([[0, 0], [1000, 1000]])],
        data={"crsid": [1, 1, 1, 1], "order": [1, 2, 3, 4], "z": [10, 0, 1, 11]}, crs=28992)
    mixed.index = mixed.index.astype(str)

    reduce_gdf_precision(linestrings)
    reduce_gdf_precision(points)
    with pytest.raises(NotImplementedError) as e:
        reduce_gdf_precision(mixed)
    assert "Mixed geometry types are not supported." in str(e.value)
