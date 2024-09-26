from os.path import abspath, dirname, join
import geopandas as gpd
import pandas as pd
import hydromt_delft3dfm.workflows.crosssections as xsec # needed to access private functions
from shapely import Point, LineString

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_set_branch_crosssections(branches):
    # xsec.set_branch_crosssections()
    pass

def test_set_xyz_crosssections():
    # test on a single branch with 2 cross-sections (id = 1 .. 2)
    branches = gpd.GeoDataFrame(data={"frictionid": ["Manning_0.023"], "frictiontype": ["Manning"], "frictionvalue": [0.023]}, geometry=[LineString([[0, 0], [1000, 1000]])])
    xyz_crosssections_1 = gpd.GeoDataFrame(
        geometry=[Point(-10, 1), Point(0, 1), Point(10, 2), Point(20, 2)],
        data={"crsid": [1, 1, 1, 1], "order": [1, 2, 3, 4], "z": [10, 0, 1, 11]})
    xyz_crosssections_2 = gpd.GeoDataFrame(
        geometry=[Point(890, 900), Point(900, 895), Point(910, 890), Point(920, 890)],
        data={"crsid": [2, 2, 2, 2], "order": [1, 2, 3, 4], "z": [9, -1, 0, 10]})
    xsec.set_xyz_crosssections(branches=branches, crosssections=pd.concat([xyz_crosssections_1, xyz_crosssections_2]))

def test_set_point_crosssections(branches, point_crosssections):
    # xsec.set_point_crosssections
    pass

def test_set_circle_crs(circle_crosssections):
    # xsec._set_circle_crs()
    pass

def test_set_rectangle_crs(rectangle_crosssections):
    # xsec._set_rectangle_crs()
    pass

def test_set_trapezoid_crs(trapezoid_crosssections):
    # xsec._set_trapezoid_crs()
    pass

def test_set_zw_crs(zw_crosssections):
    # xsec._set_zw_crs()
    pass

def test_set_yz_crs(yz_crosssections):
    # xsec._set_yz_crs()
    pass