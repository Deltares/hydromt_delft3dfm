"""Test for HydroMT plugin model class DFlowFMModel"""

import logging
from os.path import abspath, dirname, join
import pytest

from hydromt.readers import read_workflow_yaml
from hydromt_delft3dfm import DFlowFMModel

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "piave": {
        "ini": join(EXAMPLEDIR, "dflowfm_build_piave.yml"),
        "data": "artifact_data",
    },
    "local": {
        "ini": join(EXAMPLEDIR, "dflowfm_build_local.yml"),
        "data": join(EXAMPLEDIR, "data", "data_catalog_local.yaml"),
    },
}


@pytest.mark.timeout(300)  # max 5 min
@pytest.mark.slow
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
    logger.setLevel(logging.DEBUG)
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

    # somehow a successful check still has {"geoms":{}), so this should be popped
    if len(errors["geoms"]) == 0:
        errors.pop("geoms", None)
        if len(errors) == 0:
            equal = True

    assert equal, errors


@pytest.mark.timeout(120)  # max 2 min
@pytest.mark.slow
def test_model_update(tmp_path):
    # Build method options
    config = join(EXAMPLEDIR, "dflowfm_update_mesh2d_refine.yml")
    _, _, steps = read_workflow_yaml(config)

    # test update method, compare network size before and after refinement
    root = join(EXAMPLEDIR, f"dflowfm_piave")
    logger = logging.getLogger("hydromt")
    logger.setLevel(logging.DEBUG)
    model = DFlowFMModel(
        root=root,
        mode="r+",
    )
    model.read()
    netw1 = model.mesh.data.copy()
    assert len(netw1.mesh2d_node_x) == 504
    assert len(netw1.link1d2d_id) == 1734

    # set new outputdir and write mode
    outputdir = tmp_path / "dflowfm_mesh2d_refine"
    model.root.set(path=outputdir, mode="w")

    # call update() only after setting mode=w
    model.update(steps=steps)
    # alternatively call the steps from the workflow yml in Python code
    # then the data_catalog is not automatically exported on write, so do this manually
    # polygon_fn = join(EXAMPLEDIR, "data/refine.geojson")
    # model.setup_mesh2d_refine(polygon_fn=polygon_fn, steps=2)
    # model.setup_link1d2d(link_direction="1d_to_2d")
    # model.data_catalog.export_data(outputdir)

    netw2 = model.mesh.data.copy()
    assert len(netw2.mesh2d_node_x) == 767
    assert len(netw2.link1d2d_id) == 1710

    model.write()
