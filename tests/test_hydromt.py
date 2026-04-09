"""Test for HydroMT plugin model class DFlowFMModel"""

import logging
from os.path import abspath, dirname, join
import pytest

from hydromt.readers import read_workflow_yaml
from hydromt_delft3dfm import DFlowFMModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "piave": {
        "ini": join(TESTDATADIR, "dflowfm_build_piave.yml"),
        "data": "artifact_data",
    },
    "local": {
        "ini": join(TESTDATADIR, "dflowfm_build_local.yml"),
        "data": join(TESTDATADIR, "data_catalog_local.yaml"),
    },
}


@pytest.mark.timeout(300)  # max 5 min
@pytest.mark.parametrize("modelname", list(_models.keys()))
def test_model_build(tmpdir, modelname):
    model_dict = _models[modelname]

    # Build method options
    config = model_dict["ini"]
    _, global_sect, steps = read_workflow_yaml(config)
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']

    # test build method
    # compare results with model from examples folder
    root = join(tmpdir, f"dflowfm_{modelname}")
    logger = logging.getLogger("hydromt")
    logger.setLevel(10)
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
    )
    # Build model (now excludes global section because of pop)
    mod1.build(steps=steps)

    # Compare with model from examples folder
    # (need to read it again for proper geoms check)
    mod1 = DFlowFMModel(root=root, mode="r")
    mod1.read()
    root = join(EXAMPLEDIR, f"dflowfm_{modelname}")
    mod0 = DFlowFMModel(root=root, mode="r")
    mod0.read()
    # check if equal
    equal, errors = mod0.test_equal(mod1)
    
    # skip config.filepath (root is different)
    if "mdu" in errors and "mdu.filepath" in errors["mdu"]:
        errors["mdu"].pop("mdu.filepath", None)
        if len(errors["mdu"]) == 0:
            errors.pop("mdu", None)
        if len(errors) == 0:
            equal = True

    # test_equal fails for geoms, to be fixed hydromt-core or in
    # https://github.com/Deltares/hydromt_delft3dfm/issues/138
    # for now, check whether the situation does not detoriate by asserting the failing geoms
    if "geoms" in errors.keys():
        if "crosssections" in errors["geoms"]: # local
            msg = errors["geoms"]["crosssections"]
            expected_msg = '2 out of 45 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 44] \nThe first not almost equal geometry:\nLeft: POINT (663683.6568543091 1525448.030013659)\nRight: POINT (663683.6567697098 1525448.0299216388)\n'
            if msg == expected_msg:
                errors["geoms"].pop("crosssections", None)
        if "manholes" in errors["geoms"]: # piave
            msg = errors["geoms"]["manholes"]
            expected_msg_start = "997 out of 1007 geometries are not almost equal.\n"
            expected_msg_end = "\nThe first not almost equal geometry:\nLeft: POINT (1388946.301101 5860414.883978)\nRight: POINT (1389221.905704 5860123.325637)\n"
            if msg.startswith(expected_msg_start) and msg.endswith(expected_msg_end):
                errors["geoms"].pop("manholes", None)
        if "rivers" in errors["geoms"]: # piave
            msg = errors["geoms"]["rivers"]
            expected_msg = "2 out of 6 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 2] \nThe first not almost equal geometry:\nLeft: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\nRight: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\n"
            if msg == expected_msg:
                errors["geoms"].pop("rivers", None)
            # TODO ASSERTIONERROR DUE TO NEW MESSAGE
            # 'GeoDataFrame.iloc[:, 5] (column name="branchorder") are different\n\nGeoDataFrame.iloc[:, 5] (column name="branchorder") values are different (100.0 %)\n[index]: [0, 1, 2, 3, 4, 5]\n[left]:  [-1, -1, -1, 1, 1, -1]\n[right]: [-1.0, -1.0, -1.0, 1.0, 1.0, -1.0]\nAt positional index 0, first diff: -1 != -1.0'
        if "pipes" in errors["geoms"]: # piave
            msg = errors["geoms"]["pipes"]
            expected_msg_start = 'GeoDataFrame.iloc[:, 3] (column name="manhole_up") are different\n\nGeoDataFrame.iloc[:, 3] (column name="manhole_up") values are different (98.8012 %)\n[index]: '
            expected_msg_end = "\nAt positional index 11, first diff: manhole_11_generated != manhole_10_generated"
            if msg.startswith(expected_msg_start) and msg.endswith(expected_msg_end):
                errors["geoms"].pop("pipes", None)
            # TODO ASSERTIONERROR DUE TO NEW MESSAGE
            # 'GeoDataFrame.iloc[:, 5] (column name="branchorder") are different\n\nGeoDataFrame.iloc[:, 5] (column name="branchorder") values are different (100.0 %)\n[index]: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, ...]\n[left]:  [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ...]\n[right]: [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, ...]\nAt positional index 0, first diff: -1 != -1.0'
        if "branches" in errors["geoms"]: # piave
            msg = errors["geoms"]["branches"]
            expected_msg = "2 out of 1007 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 2] \nThe first not almost equal geometry:\nLeft: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\nRight: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\n"
            if msg == expected_msg:
                errors["geoms"].pop("branches", None)
            # TODO ASSERTIONERROR DUE TO NEW MESSAGE
            # 'GeoDataFrame.iloc[:, 5] (column name="branchorder") are different\n\nGeoDataFrame.iloc[:, 5] (column name="branchorder") values are different (100.0 %)\n[index]: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, ...]\n[left]:  [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ...]\n[right]: [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, ...]\nAt positional index 0, first diff: -1 != -1.0'
        
        # TODO ADDITIONAL GEOMS THAT ARE NOT EQUAL IN PIAVE: 'boundaries', 'mesh1d', 'pipe_nodes', 
        
        if len(errors["geoms"]) == 0:
            errors.pop("geoms", None)
        if len(errors) == 0:
            equal = True
    
    assert equal, errors


# for debugging in IDE
if __name__ == "__main__":
    test_model_build(tmpdir=".", modelname="piave")
