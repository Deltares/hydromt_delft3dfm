from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
import numpy as np
from pathlib import Path

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")
TOLERANCE = 1e-6


def test_write_empty_model(tmpdir):
    """
    failed before: https://github.com/Deltares/hydromt_delft3dfm/issues/246
    """
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mod1.write()


def test_write_readonlymode(tmpdir, caplog):
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mod1.write()
    
    mod2 = DFlowFMModel(
        root=root,
        mode="r",
        crs=3857,
    )
    mod2.write()
    assert "Cannot write in read-only mode" in caplog.text


def test_read_write_config_empty_paths(tmpdir):
    # Instantiate an empty model
    root = join(tmpdir, "dflowfm_example")
    model1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    # Get the mdu settings
    model1.mdu.read()
    # Check whether the path is an emtpy string
    assert model1.mdu.data["output"]["outputdir"] == ""
    
    # write the model to read it again
    model1.write()
    model2 = DFlowFMModel(root=root, mode="r", crs=3857)
    # Get the mdu settings
    model2.mdu.read()
    # Check whether the path is an emtpy string
    # TODO: should be an empty string: https://github.com/Deltares/HYDROLIB-core/issues/703
    # then update this test: https://github.com/Deltares/hydromt_delft3dfm/issues/148
    # and update the reference mdu files for piave and local
    assert model2.mdu.data["output"]["outputdir"] == Path(".")


def test_setup_mesh2d_refine(tmpdir):
    # get dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), crs=3857, mode="r")
    mesh2d = model.mesh.get_mesh('mesh2d')
    assert mesh2d.face_coordinates.shape == (460, 2)
    assert mesh2d.edge_coordinates.shape == (963, 2)
    mesh1d = model.mesh.get_mesh('mesh1d')
    assert mesh1d.edge_coordinates.shape == (1732, 2)

    # refine and assert
    model.setup_mesh2d_refine(polygon_fn=join(EXAMPLEDIR, "data","refine.geojson"))
    mesh2d = model.mesh.get_mesh('mesh2d')
    assert mesh2d.face_coordinates.shape == (656, 2)
    assert mesh2d.edge_coordinates.shape == (1306, 2)
    mesh1d = model.mesh.get_mesh('mesh1d')
    assert mesh1d.edge_coordinates.shape == (1732, 2)


def test_setup_channels(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")

    # setup_channels
    region = {'geom': join(TESTDATADIR, "local_data","1D_extent.geojson")}
    channels_fn = join(TESTDATADIR, "local_data","1D_rivers.geojson")
    crosssections_fn = join(TESTDATADIR, "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region,
        channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )


def test_setup_retentions(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    
    # setup_retentions
    retentions_fn = join(TESTDATADIR, "local_data","retention_ponds.geojson")
    # Add 1 retention pond, should be included with snap_offset = 200
    model.setup_retentions(retentions_fn=retentions_fn, snap_offset=200)
    assert len(model.geoms.data["retentions"]) == 1


def test_setup_bridges(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    
    # first add channels to obtain friction values for branches
    # see also https://github.com/Deltares/hydromt_delft3dfm/issues/168
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
    assert len(model.geoms.data['bridges']) == 2 


def test_setup_culverts(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")

    # first add channels to obtain friction values for branches
    # see also https://github.com/Deltares/hydromt_delft3dfm/issues/168
    region = {'geom': join(TESTDATADIR, "local_data","1D_extent.geojson")}
    channels_fn = join(TESTDATADIR, "local_data","1D_rivers.geojson")
    crosssections_fn = join(TESTDATADIR, "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )

    # setup culverts (total of 1 culvert)
    culverts_fn = join(TESTDATADIR, "local_data","culverts.geojson")
    model.setup_culverts(culverts_fn=culverts_fn)
    assert len(model.geoms.data['culverts']) == 1


def test_write_structures(tmpdir):
    """
    failed before for dflowfm_local model due to nan values in gdf
    https://github.com/Deltares/hydromt_delft3dfm/issues/150
    """
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    
    # indirectly call hidden write_structures() method
    model.geoms.write(write_mesh_gdf=False)


def test_inifield_add_raster_data_from_rasterdataset(tmpdir):
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), crs=3857, mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    raster_fn = join(TESTDATADIR, "local_data","frictioncoefficient.tif")
    variable = 'roughness_manning'
    variables = [variable]
    model.inifield.add_raster_data_from_rasterdataset(raster_fn, variables)

    roughness_values = np.unique(model.inifield.data[variable]).tolist()
    expected_values = [-999.0, 0.025, 0.044, 0.050, 0.055]
    assert np.allclose(roughness_values, expected_values, atol=TOLERANCE)