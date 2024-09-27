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


def test_setup_links1d2d_add_links(tmpdir):
    # Instantiate an empty model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.set_root(tmpdir, mode="w")
    
    # TODO: this test takes >6 seconds, speed it up by making it more unittest-like

    # TODO: add checks with assertions, but at the moment nothing much seems to change
    # compared to the initial model (eg the shape of 1d2dlinks is still (102,2))

    # these lines below are still useful to increase code coverage
    # at least the functions are executed and raise no errors
    model.setup_link1d2d(link_direction="1d_to_2d", link_type="embedded")
    model.setup_link1d2d(link_direction="2d_to_1d", link_type="embedded")
    model.setup_link1d2d(link_direction="1d_to_2d", link_type="lateral")
    # TODO: below does not work with dist_factor enabled
    model.setup_link1d2d(link_direction="2d_to_1d", link_type="lateral", dist_factor="None")
    
    # much of the validation is done upon writing, so maybe write the result also
    model.write_mesh()
    model.write()
