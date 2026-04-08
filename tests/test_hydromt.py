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
    # check if equal, geoms are temporarily skipped
    # https://github.com/Deltares/hydromt_delft3dfm/issues/138
    equal, errors = mod0._test_equal(mod1, skip_component=["geoms"])
    # skip config.filepath (root is different)
    if "config.filepath" in errors:
        errors.pop("config.filepath", None)
        if len(errors) == 0:
            equal = True
    assert equal, errors


def test_model_build_local_code(tmp_path):
    """
    A python code version of the local model, to make debugging easier.
    This test can be removed once the entire hydromt_delft3dfm is properly
    covered by unittests.
    """
    model_dict = _models["local"]
    # Build method options
    config = model_dict["ini"]
    _, global_sect, steps = read_workflow_yaml(config)
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']
    # initialize model
    logger = logging.getLogger("hydromt")
    logger.setLevel(10)
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
    )
    # build model via steps corresponding to yml order
    for step_dict in steps:
        step, kwargs = next(iter(step_dict.items()))
        method = getattr(model, step)
        print(f"CALLING model.{step} with input:", kwargs)
        method(**kwargs)


def test_model_build_piave_code(tmp_path):
    """
    A python code version of the piave model, to make debugging easier.
    This test can be removed once the entire hydromt_delft3dfm is properly
    covered by unittests.
    """
    model_dict = _models["piave"]
    # Build method options
    config = model_dict["ini"]
    _, global_sect, steps = read_workflow_yaml(config)
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']
    # initialize model
    logger = logging.getLogger("hydromt")
    logger.setLevel(10)
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
    )
    # build model via steps corresponding to yml order
    for step_dict in steps:
        step, kwargs = next(iter(step_dict.items()))
        method = getattr(model, step)
        print(f"CALLING model.{step} with input:", kwargs)
        method(**kwargs)
