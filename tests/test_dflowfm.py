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


def test_setup_channels(tmpdir):
    """
    to increase code coverage, we are not actually doing anything here since all options fail
    nevertheless a convenient setup for a sort of unit test, 
    so good to keep somewhere and improve in the future
    """
    # Instantiate an empty model
    model = DFlowFMModel(
        root=join(EXAMPLEDIR, "dflowfm_piave"), 
        mode="r", 
        data_libs=[join(TESTDATADIR, "test_data.yaml")]
    )
    model.read()
    model.set_root(tmpdir, mode="w")
    region = {'geom': join(TESTDATADIR, "local_data","1D_extent.geojson")}
    channels_fn = join(TESTDATADIR, "local_data","1D_rivers.geojson")
    crosssections_fn = join(TESTDATADIR, "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )


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
