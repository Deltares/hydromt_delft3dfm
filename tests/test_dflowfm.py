import pytest
from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
from pathlib import Path

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_read_write_config_empty_paths(tmpdir):
    # Instantiate an empty model
    dir_root = join(EXAMPLEDIR, "dflowfm_piave")
    dir_model = join(tmpdir, "dflowfm_piave")
    import shutil
    shutil.copytree(dir_root, dir_model)
    model = DFlowFMModel(root=dir_model, mode="r+")
    # Get the mdu settings
    model.read_config()
    # Check whether the path is an emtpy string
    # TODO: we temporarly put . in the example mdu, so this is now also here
    assert model.config["output"]["outputdir"] == Path(".") 
    
    # write the mdu to read again
    model.write_config()
    # re-read the model
    model2 = DFlowFMModel(root=dir_model, mode="r")
    # Get the mdu settings
    model2.read_config()
    # Check whether the path is an emtpy string
    # TODO: should be an empty string: https://github.com/Deltares/HYDROLIB-core/issues/703
    # then update this test: https://github.com/Deltares/hydromt_delft3dfm/issues/148
    assert model2.config["output"]["outputdir"] == Path(".")


def test_setup_mesh2d_refine(tmpdir):
    # get dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r")
    mesh2d = model.get_mesh('mesh2d')
    assert mesh2d.face_coordinates.shape == (460, 2)
    assert mesh2d.edge_coordinates.shape == (963, 2)
    mesh1d = model.get_mesh('mesh1d')
    assert mesh1d.edge_coordinates.shape == (1732, 2)

    # refine and assert
    model.setup_mesh2d_refine(polygon_fn=join(EXAMPLEDIR, "data","refine.geojson"))
    mesh2d = model.get_mesh('mesh2d')
    assert mesh2d.face_coordinates.shape == (656, 2)
    assert mesh2d.edge_coordinates.shape == (1306, 2)
    mesh1d = model.get_mesh('mesh1d')
    assert mesh1d.edge_coordinates.shape == (1732, 2)


def test_setup_channels(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(
        root=join(EXAMPLEDIR, "dflowfm_local"), 
        mode="r", 
        data_libs=[join(TESTDATADIR, "test_data.yaml")]
    )
    model.read()
    model.set_root(tmpdir, mode="w")

    # setup_channels
    region = {'geom': join(TESTDATADIR, "local_data","1D_extent.geojson")}
    channels_fn = join(TESTDATADIR, "local_data","1D_rivers.geojson")
    crosssections_fn = join(TESTDATADIR, "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )


def test_setup_bridges(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(
        root=join(EXAMPLEDIR, "dflowfm_local"), 
        mode="r", 
        data_libs=[join(TESTDATADIR, "test_data.yaml")]
    )
    model.read()
    model.set_root(tmpdir, mode="w")
    
    # first add channels to obtain friction values for branches
    region = {'geom': join(TESTDATADIR, "local_data","1D_extent.geojson")}
    channels_fn = join(TESTDATADIR, "local_data","1D_rivers.geojson")
    crosssections_fn = join(TESTDATADIR, "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )

    # setup bridges (total of 2 bridges)
    bridges_fn = join(TESTDATADIR, "local_data","bridges.geojson")
    model.setup_bridges(bridges_fn=bridges_fn)
    assert len(model.geoms['bridges']) == 2 

def test_setup_culverts(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(
        root=join(EXAMPLEDIR, "dflowfm_local"), 
        mode="r", 
        data_libs=[join(TESTDATADIR, "test_data.yaml")]
    )
    model.read()
    model.set_root(tmpdir, mode="w")

    # setup culverts
    culverts_fn = join(TESTDATADIR, "local_data","culverts.geojson")
    model.setup_culverts(culverts_fn=culverts_fn)


def test_write_structures(tmpdir):
    """
    failed before for dflowfm_local model due to nan values in gdf
    https://github.com/Deltares/hydromt_delft3dfm/issues/150
    """
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.set_root(tmpdir, mode="w")
    
    # indirectly call hidden write_structures() method
    model.write_geoms(write_mesh_gdf=False)
