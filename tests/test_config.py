from hydromt_delft3dfm import DFlowFMModel

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

def test_model_config(tmpdir):
    # Instantiate an empty model
    model = DFlowFMModel(root=tmpdir, mode="w")
    # Get the template mdu
    model.read_config()
    # Check an option
    assert model.config["output"]["outputdir"] = "output"
    model.write_config()

    model2 = DFlowFMModel(root=tmpdir, mode="r+")
    model2.read_config()

def test_read_config(tmpdir):
    # Instantiate an empty model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r+")
    # Get the template mdu
    model.read_config()
    # Check an option
    assert model.config["output"]["outputdir"] = "output"

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
