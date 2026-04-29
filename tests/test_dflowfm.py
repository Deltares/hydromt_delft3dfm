from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
import numpy as np
from pathlib import Path
import pytest

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


def test_init_dflowfmmodel_mode_write_crs_none(tmpdir):
    """
    tests whether the crs is parsed properly
    https://github.com/Deltares/hydromt_delft3dfm/issues/247
    """
    root = join(tmpdir, "dflowfm_example")
    with pytest.raises(ValueError) as e:
        _ = DFlowFMModel(root=root, mode="w")
    assert "crs argument cannot be None with mode" in str(e.value)
    with pytest.raises(ValueError) as e:
        _ = DFlowFMModel(root=root, mode="w+")
    assert "crs argument cannot be None with mode" in str(e.value)


def test_init_dflowfmmodel_mode_read_crs_none():
    """
    tests whether the crs is parsed properly
    https://github.com/Deltares/hydromt_delft3dfm/issues/247
    """
    root = join(EXAMPLEDIR, "dflowfm_local")
    model1 = DFlowFMModel(root=root, mode="r")
    model2 = DFlowFMModel(root=root, mode="r+")
    assert model1.crs.to_epsg() == 32647
    assert model2.crs.to_epsg() == 32647


def test_init_dflowfmmodel_mode_read_crs_notnone(tmpdir):
    """
    tests whether the crs is parsed properly
    https://github.com/Deltares/hydromt_delft3dfm/issues/247
    """
    root = join(EXAMPLEDIR, "dflowfm_local")
    model1 = DFlowFMModel(root=root, mode="r", crs=None)
    model2 = DFlowFMModel(root=root, mode="r", crs=None)
    assert model1.crs.to_epsg() == 32647
    assert model2.crs.to_epsg() == 32647

    # TODO: the model crs is actually 32647, but is overwritten here
    # https://github.com/Deltares/hydromt_delft3dfm/issues/119
    model3 = DFlowFMModel(root=root, mode="r", crs=4326)
    model4 = DFlowFMModel(root=root, mode="r", crs=4326)
    assert model3.crs.to_epsg() == 4326
    assert model4.crs.to_epsg() == 4326


def test_write_readonlymode(tmpdir, caplog):
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    # TODO: cannot read the crs of a model if no mesh was included
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    mod1.write()
    
    mod2 = DFlowFMModel(root=root, mode="r")
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

    # TODO: cannot read the crs of a model if no mesh was included
    model1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )

    # write the model to read it again
    model1.write()
    model2 = DFlowFMModel(root=root, mode="r")
    # Get the mdu settings
    model2.mdu.read()
    # Check whether the path is an emtpy string
    # TODO: should be an empty string: https://github.com/Deltares/HYDROLIB-core/issues/703
    # then update this test: https://github.com/Deltares/hydromt_delft3dfm/issues/148
    # and update the reference mdu files for piave and local
    assert model2.mdu.data["output"]["outputdir"] == Path(".")


def test_setup_mesh2d_refine(tmpdir):
    # get dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_piave"), mode="r")
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
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
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
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    
    # setup_retentions
    retentions_fn = join(TESTDATADIR, "local_data","retention_ponds.geojson")
    # Add 1 retention pond, should be included with snap_offset = 200
    model.setup_retentions(retentions_fn=retentions_fn, snap_offset=200)
    assert len(model.geoms.data["retentions"]) == 1


def test_setup_bridges(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
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
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
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
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    
    # indirectly call hidden write_structures() method
    model.geoms.write(write_mesh_gdf=False)


def test_inifield_add_raster_data_from_rasterdataset(tmpdir):
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")
    raster_fn = join(TESTDATADIR, "local_data","frictioncoefficient.tif")
    variable = 'roughness_manning'
    variables = [variable]
    model.inifield.add_raster_data_from_rasterdataset(raster_fn, variables)

    roughness_values = np.unique(model.inifield.data[variable]).tolist()
    expected_values = [-999.0, 0.025, 0.044, 0.050, 0.055]
    assert np.allclose(roughness_values, expected_values, atol=TOLERANCE)