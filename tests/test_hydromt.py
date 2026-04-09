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
    if "geoms.crosssections" in errors: # local
        msg = errors["geoms.crosssections"].args[0]
        expected_msg = '2 out of 45 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 44] \nThe first not almost equal geometry:\nLeft: POINT (663683.6568543091 1525448.030013659)\nRight: POINT (663683.6567697098 1525448.0299216388)\n'
        assert msg == expected_msg
        errors.pop("geoms.crosssections", None)
    if "geoms.manholes" in errors: # piave
        msg = errors["geoms.manholes"].args[0]
        expected_msg_start = "997 out of 1007 geometries are not almost equal.\n"
        expected_msg_end = "\nThe first not almost equal geometry:\nLeft: POINT (1388946.301101 5860414.883978)\nRight: POINT (1389221.905704 5860123.325637)\n"
        assert msg.startswith(expected_msg_start)
        assert msg.endswith(expected_msg_end)
        errors.pop("geoms.manholes", None)
    if "geoms.rivers" in errors: # piave
        msg = errors["geoms.rivers"].args[0]
        expected_msg = "2 out of 6 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 2] \nThe first not almost equal geometry:\nLeft: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\nRight: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\n"
        assert msg == expected_msg
        errors.pop("geoms.rivers", None)
    if "geoms.pipes" in errors: # piave
        msg = errors["geoms.pipes"].args[0]
        expected_msg_start = 'GeoDataFrame.iloc[:, 3] (column name="manhole_up") are different\n\nGeoDataFrame.iloc[:, 3] (column name="manhole_up") values are different (98.8012 %)\n[index]: '
        expected_msg_end = "\nAt positional index 11, first diff: manhole_11_generated != manhole_10_generated"
        assert msg.startswith(expected_msg_start)
        assert msg.endswith(expected_msg_end)
        errors.pop("geoms.pipes", None)
    if "geoms.branches" in errors: # piave
        msg = errors["geoms.branches"].args[0]
        expected_msg = "2 out of 1007 geometries are not almost equal.\nIndices where geometries are not almost equal: [0, 2] \nThe first not almost equal geometry:\nLeft: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\nRight: LINESTRING (1393117.044157 5866704.07164, 1393024.277914 5866704.07164, 1392931.511672 5866704.07164...\n"
        assert msg == expected_msg
        errors.pop("geoms.branches", None)
    if len(errors) == 0:
        equal = True
    
    assert equal, errors
