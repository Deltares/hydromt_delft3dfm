from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
from hydromt_delft3dfm.mesh_utils import hydrolib_network_from_mesh

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_hydrolib_network_from_mesh(tmpdir):
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r")
    model.read()
    model.set_root(tmpdir, mode="w")

    # check if all expected grids and link1d2 are present in mesh
    assert set(model.mesh_names) == set(['mesh1d', 'network1d', 'mesh2d'])
    assert "link1d2d" in model.mesh.data_vars
    assert model.mesh["link1d2d"].shape == (1734,2)

    # check 1d2dlinks are properly converted to hydrolib-core network/link1d2d instance
    network = hydrolib_network_from_mesh(model.mesh)
    link1d2d_arr = network._link1d2d.link1d2d
    assert link1d2d_arr.shape == (1734, 2)

    # call this since it does some checks
    model.write_mesh()


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

    # TODO: add checks with assertions, but at the moment nothing much seems to change compared to the initial model
    # at least not with respect to 1d2dlinks shape

    # these lines below are still useful to increase code coverage
    # at least the functions are executed and raise no errors
    model.setup_link1d2d(link_direction="1d_to_2d", link_type="embedded")
    model.setup_link1d2d(link_direction="2d_to_1d", link_type="embedded")
    model.setup_link1d2d(link_direction="1d_to_2d", link_type="lateral")
    # TODO: uncomment this last case, now raises `TypeError: unsupported operand type(s) for -: 'float' and 'Point'` 
    # (from workflows/mesh.py line 688)
    # model.setup_link1d2d(link_direction="2d_to_1d", link_type="lateral")
    

    # much of the validation is done upon writing, so maybe write the result also
    model.write_mesh()
    model.write()
