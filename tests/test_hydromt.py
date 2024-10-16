"""Test for hydromt plugin model class DFlowFMModel"""

import pdb
from os.path import abspath, dirname, join

import pytest
from hydromt.cli.cli_utils import parse_config
from hydromt.log import setuplog

from hydromt_delft3dfm import DFlowFMModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "piave": {
        "example": "dflowfm_piave",
        "ini": "dflowfm_build.yml",
        "data": "artifact_data",
        "snap_offset": 25,
        "crs": 3857,  # global section needs to be passed to build as arguments
    },
    "local": {
        "example": "dflowfm_local",
        "ini": "dflowfm_build_local.yml",
        "data": join(TESTDATADIR, "test_data.yaml"),
        "snap_offset": 25,
        "crs": 32647,  # global section needs to be passed to build as arguments
    },
}


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_class(model):
    _model = _models[model]
    # read model in examples folder
    root = join(EXAMPLEDIR, _model["example"])
    mod = DFlowFMModel(root=root, mode="r")
    mod.read()
    # run test_model_api() method
    non_compliant_list = mod._test_model_api()
    assert len(non_compliant_list) == 0


@pytest.mark.timeout(300)  # max 5 min
@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_build(tmpdir, model):
    _model = _models[model]
    # test build method
    # compare results with model from examples folder
    root = join(tmpdir, model)
    logger = setuplog(__name__, join(root, "hydromt.log"), log_level=10)
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        data_libs=[_model["data"]],
        network_snap_offset=_model["snap_offset"],
        crs=_model["crs"],
        openwater_computation_node_distance=40,
        logger=logger,
    )
    # Build method options
    config = join(TESTDATADIR, _model["ini"])
    opt = parse_config(config)
    # pop global section and assert values with init
    global_sect = opt.pop('global')
    assert mod1._crs == global_sect['crs']
    assert mod1._network_snap_offset == global_sect['network_snap_offset']
    assert (mod1._openwater_computation_node_distance == 
            global_sect['openwater_computation_node_distance'])
    # Build model
    mod1.build(opt=opt)
    # Check if model is api compliant
    non_compliant_list = mod1._test_model_api()
    assert len(non_compliant_list) == 0

    # Compare with model from examples folder
    # (need to read it again for proper geoms check)
    mod1 = DFlowFMModel(root=root, mode="r", logger=logger)
    mod1.read()
    root = join(EXAMPLEDIR, _model["example"])
    mod0 = DFlowFMModel(root=root, mode="r")
    mod0.read()
    # check if equal
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
    logger = setuplog(__name__, join(tmp_path, "hydromt.log"), log_level=10)
    _model = _models["local"]
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[_model["data"]],
        network_snap_offset=_model["snap_offset"],
        crs=_model["crs"],
        openwater_computation_node_distance=40,
        logger=logger
    )

    config = join(TESTDATADIR, _model["ini"])
    opt = parse_config(config)
    model.setup_rivers(**opt['setup_rivers'])
    model.setup_rivers(**opt['setup_rivers1'])
    model.setup_pipes(**opt['setup_pipes'])
    model.setup_manholes(**opt['setup_manholes'])
    model.setup_bridges(**opt['setup_bridges'])
    model.setup_culverts(**opt['setup_culverts'])
    model.setup_1dboundary(**opt['setup_1dboundary'])
    model.setup_1dlateral_from_points(**opt['setup_1dlateral_from_points'])
    model.setup_1dlateral_from_polygons(**opt['setup_1dlateral_from_polygons'])
    model.setup_mesh2d(**opt['setup_mesh2d'])
    # model.setup_mesh2d_refine(**opt['setup_mesh2d_refine'])
    model.setup_maps_from_rasterdataset(**opt['setup_maps_from_rasterdataset'])
    model.setup_2dboundary(**opt['setup_2dboundary'])
    model.setup_rainfall_from_uniform_timeseries(**opt['setup_rainfall_from_uniform_timeseries'])
    model.setup_link1d2d(**opt['setup_link1d2d'])


def test_model_build_piave_code(tmp_path):
    """
    A python code version of the piave model, to make debugging easier.
    This test can be removed once the entire hydromt_delft3dfm is properly
    covered by unittests.
    """
    logger = setuplog(__name__, join(tmp_path, "hydromt.log"), log_level=10)
    _model = _models["piave"]
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[_model["data"]],
        network_snap_offset=_model["snap_offset"],
        crs=_model["crs"],
        openwater_computation_node_distance=40,
        logger=logger
    )

    config = join(TESTDATADIR, _model["ini"])
    opt = parse_config(config)
    model.setup_rivers_from_dem(**opt['setup_rivers_from_dem'])
    model.setup_pipes(**opt['setup_pipes'])
    model.setup_manholes(**opt['setup_manholes'])
    model.setup_1dboundary(**opt['setup_1dboundary'])
    model.setup_mesh2d(**opt['setup_mesh2d'])
    model.setup_maps_from_rasterdataset(**opt['setup_maps_from_rasterdataset'])
    model.setup_maps_from_raster_reclass(**opt['setup_maps_from_raster_reclass'])
    model.setup_link1d2d(**opt['setup_link1d2d'])
