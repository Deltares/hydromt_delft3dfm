from os.path import abspath, dirname, join
import geopandas as gpd
import pandas as pd
import hydromt_delft3dfm.workflows.crosssections as xsec # needed to access private functions
from shapely import Point, LineString

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_set_branch_crosssections():
    # xsec.set_branch_crosssections(branches)
    pass

def test_set_xyz_crosssections():
    # test on a single branch with 3 cross-sections (id = 1 .. 3)
    branches = gpd.GeoDataFrame(data={"frictionid": ["Manning_0.023"], "frictiontype": ["Manning"], "frictionvalue": [0.023]}, geometry=[LineString([[0, 0], [1000, 1000]])], crs=28992)
    
    xyz_crosssections_1 = gpd.GeoDataFrame(
        geometry=[Point(-10, 1), Point(0, 1), Point(10, 2), Point(20, 2)],
        data={"crsid": [1, 1, 1, 1], "order": [1, 2, 3, 4], "z": [10, 0, 1, 11]}, crs=28992)
    xyz_crosssections_2 = gpd.GeoDataFrame(
        geometry=[Point(890, 900), Point(900, 895), Point(910, 890), Point(920, 890)],
        data={"crsid": [2, 2, 2, 2], "order": [1, 2, 3, 4], "z": [9, -1, 0, 10]}, crs=28992)
    # cross-section 3 deliberately has too few points
    xyz_crosssections_3 = gpd.GeoDataFrame(
        geometry=[Point(995, 995), Point(1005, 1005)],
        data={"crsid": [3, 3], "order": [1, 2], "z": [8, 8]}, crs=28992)
    
    xyz_crosssections = pd.concat([xyz_crosssections_1, xyz_crosssections_2, xyz_crosssections_3]).reset_index(drop=True)
    crosssections = xsec.set_xyz_crosssections(branches=branches, crosssections=xyz_crosssections)

    assert len(crosssections) == 2 # check if the cross-section containing only 2 points is dropped
    #TODO more checks can be added 

def test_set_point_crosssections():
    # xsec.set_point_crosssections(branches, point_crosssections)
    pass

def test_set_circle_crs():
    # xsec._set_circle_crs(circle_crosssections)
    pass

def test_set_rectangle_crs():
    # xsec._set_rectangle_crs(rectangle_crosssections)
    pass

def test_set_trapezoid_crs():
    # xsec._set_trapezoid_crs(trapezoid_crosssections)
    pass

def test_set_zw_crs():
    # xsec._set_zw_crs(zw_crosssections)
    pass

def test_set_yz_crs():
    # xsec._set_yz_crs(yz_crosssections)
    pass