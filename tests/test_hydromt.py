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
        "ini": join(TESTDATADIR, "dflowfm_build_piave.yml"),
        "data": "artifact_data",
    },
    "local": {
        "ini": join(TESTDATADIR, "dflowfm_build_local.yml"),
        "data": join(TESTDATADIR, "data_catalog_local.yaml"),
    },
}


@pytest.mark.parametrize("modelname", list(_models.keys()))
def test_model_class(modelname):
    # read model in examples folder
    root = join(EXAMPLEDIR, f"dflowfm_{modelname}")
    mod = DFlowFMModel(root=root, mode="r")
    mod.read()
    # run test_model_api() method
    non_compliant_list = mod._test_model_api()
    assert len(non_compliant_list) == 0


@pytest.mark.timeout(300)  # max 5 min
@pytest.mark.parametrize("modelname", list(_models.keys()))
def test_model_build(tmpdir, modelname):
    model_dict = _models[modelname]

    # Build method options
    config = model_dict["ini"]
    opt = parse_config(config)
    # pop global section and get model init arguments
    global_sect = opt.pop('global')
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']

    # test build method
    # compare results with model from examples folder
    root = join(tmpdir, f"dflowfm_{modelname}")
    logger = setuplog(__name__, join(root, "hydromt.log"), log_level=10)
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
        logger=logger,
    )
    # Build model (now excludes global section because of pop)
    mod1.build(opt=opt)
    # Check if model is api compliant
    non_compliant_list = mod1._test_model_api()
    assert len(non_compliant_list) == 0

    # Compare with model from examples folder
    # (need to read it again for proper geoms check)
    mod1 = DFlowFMModel(root=root, mode="r", logger=logger)
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
    opt = parse_config(config)
    # pop global section and get model init arguments
    global_sect = opt.pop('global')
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']
    # initialize model
    logger = setuplog(__name__, join(tmp_path, "hydromt.log"), log_level=10)
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
        logger=logger
    )
    # build model via steps corresponding to yml order
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
    model_dict = _models["piave"]
    # Build method options
    config = model_dict["ini"]
    opt = parse_config(config)
    # pop global section and get model init arguments
    global_sect = opt.pop('global')
    crs = global_sect['crs']
    network_snap_offset = global_sect['network_snap_offset']
    openwater_computation_node_distance = global_sect['openwater_computation_node_distance']
    # initialize model
    logger = setuplog(__name__, join(tmp_path, "hydromt.log"), log_level=10)
    model = DFlowFMModel(
        root=tmp_path,
        mode="w",
        data_libs=[model_dict["data"]],
        network_snap_offset=network_snap_offset,
        crs=crs,
        openwater_computation_node_distance=openwater_computation_node_distance,
        logger=logger
    )
    # build model via steps corresponding to yml order
    model.setup_rivers_from_dem(**opt['setup_rivers_from_dem'])
    model.setup_pipes(**opt['setup_pipes'])
    model.setup_manholes(**opt['setup_manholes'])
    model.setup_1dboundary(**opt['setup_1dboundary'])
    model.setup_mesh2d(**opt['setup_mesh2d'])
    model.setup_maps_from_rasterdataset(**opt['setup_maps_from_rasterdataset'])
    model.setup_maps_from_raster_reclass(**opt['setup_maps_from_raster_reclass'])
    model.setup_rainfall_from_constant(**opt["setup_rainfall_from_constant"])
    model.setup_link1d2d(**opt['setup_link1d2d'])
