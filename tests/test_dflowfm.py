from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
from pathlib import Path

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
