from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_setup_inks1d2d_add_links(tmpdir):
    """
    to increase code coverage, we are not actually doing anything here since all options fail
    nevertheless a convenient setup for a sort of unit test, 
    so good to keep somewhere and improve in the future
    """
    # Instantiate an empty model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r")
    model.read()
    model.set_root(tmpdir, mode="w")
    
    # TODO: make fixture of model?
    # TODO: uncomment last case
    # TODO: add checks with assertions
    
    model.setup_link1d2d(link_direction= "1d_to_2d", link_type="embedded")
    model.setup_link1d2d(link_direction= "2d_to_1d", link_type="embedded")
    model.setup_link1d2d(link_direction= "1d_to_2d", link_type="lateral")
    # model.setup_link1d2d(link_direction= "2d_to_1d", link_type="lateral")

    # much of the validation is done upon writing, so maybe write the result also
    model.write_mesh()
    model.write()
