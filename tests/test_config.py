from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


# from pathlib import Path; tmpdir = Path(r"c:\Users\veenstra\Downloads\hydromt_tmpdir")
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
    assert model.config["output"]["outputdir"] == ""
    
    # write the mdu to read again
    model.write_config()
    # re-read the model
    model2 = DFlowFMModel(root=dir_model, mode="r")
    # Get the mdu settings
    model2.read_config()
    # Check whether the path is an emtpy string
    # TODO: should be an empty string: https://github.com/Deltares/HYDROLIB-core/issues/703
    from pathlib import Path
    assert model2.config["output"]["outputdir"] == Path(".")


def test_setup_link1d2d(tmpdir):
    # Instantiate an empty model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r")
    model.read()
    model.set_root(tmpdir, mode="w")
    model.setup_link1d2d(link_direction= "1d_to_2d", link_type="embedded")
    model.setup_link1d2d(link_direction= "2d_to_1d", link_type="embedded")
    model.setup_link1d2d(link_direction= "1d_to_2d", link_type="lateral")
    # add checks with assertions
    model.write_mesh()
    model.write()
