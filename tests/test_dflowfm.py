from os.path import abspath, dirname, join
from os import makedirs, rename
from hydromt_delft3dfm import DFlowFMModel
import numpy as np
from pathlib import Path
import pytest
import shutil
import xugrid as xu

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")
TOLERANCE = 1e-6


def test_write_read_empty_model(tmpdir):
    """
    writing a model without a mesh is prohibited since
    https://github.com/Deltares/hydromt_delft3dfm/issues/270
    """
    crs = 3857
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(root=root, mode="w", crs=crs)
    with pytest.raises(RuntimeError) as e:
        mod1.write()
    expected_error = "hydromt_delft3dfm cannot write a model without a mesh/network"
    assert expected_error in str(e.value)


def test_read_empty_root_folder(tmpdir):
    """
    if writing the model fails like in test_write_read_empty_model, a root folder is
    still created (without any files). Give a proper error in this case. Added in
    https://github.com/Deltares/hydromt_delft3dfm/issues/270
    """
    # pointing to non-existent folder fails during the initialisation of the hydromt
    # parent model.
    root = join(tmpdir, "dflowfm_example")
    with pytest.raises(OSError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "model root not found at" in str(e.value)

    # create the root directory
    makedirs(root, exist_ok=False)
    # point to an empty folder fails because the network file cannot be found.
    # This only happens because crs=None, so DFlowFMModel._check_crs() activates the
    # reading of the mesh during the initialisation of the DFlowFMModel.
    with pytest.raises(ValueError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "hydromt_delft3dfm cannot read a model without a mesh/network." in str(e.value)

    # when providing a crs, DFlowFMModel._check_crs() does not activate the reading the
    # mesh during initialisation of the DFlowFMModel, but then the reading fails later
    # when calling something that tries to read the mesh.
    # TODO: it might be better to also fail on a missing mdu file immediately.
    mod3 = DFlowFMModel(root=root, mode="r", crs=4326)
    with pytest.raises(ValueError) as e:
        mod3.read()
    assert "hydromt_delft3dfm cannot read a model without a mesh/network." in str(e.value)


def test_write_read_mesh_model_different_mdu_mesh_path(tmpdir):
    """
    A model consists of at least a dimr (optional on read), a mdu and a network.
    However, the paths can be different for each model, so check if write/read
    still works when using non-default values.
    And check if the crs is preserved with write/read.
    """
    crs = 3857
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(root=root, mode="w", crs=crs, mdu_filename="folder/nonstandard.mdu")
    geom_bedlevuni = -983
    geom_netfile = "network_file_net.nc"
    mod1.setup_config(**{
        "geometry.bedlevuni": geom_bedlevuni,
        "geometry.netfile": geom_netfile,
    })
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    mod1.write()
    assert mod1.mdu.data["geometry"]["bedlevuni"] == geom_bedlevuni
    assert str(mod1.mdu.data["geometry"]["netfile"]) == geom_netfile
    # components are only available after write
    assert str(mod1.dimr.data.component[0].workingDir) == 'folder'
    assert str(mod1.dimr.data.component[0].inputFile) == 'nonstandard.mdu'

    # read in the model to see if all changes are preserved
    mod2 = DFlowFMModel(root=root, mode="r")
    assert mod1.crs.to_epsg() == crs
    assert mod2.crs.to_epsg() == crs

    # check if the updated mdu was read and not newly initialized
    assert mod2.mdu.data["geometry"]["bedlevuni"] == geom_bedlevuni
    assert str(mod2.mdu.data["geometry"]["netfile"]) == geom_netfile
    assert str(mod2.dimr.data.component[0].workingDir) == 'folder'
    assert str(mod2.dimr.data.component[0].inputFile) == 'nonstandard.mdu'


def test_write_read_model_without_dimr(tmpdir):
    crs = 3857
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(root=root, mode="w", crs=crs)
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    mod1.setup_config(**{"geometry.bedlevuni":-983})
    mod1.write()

    # remove dimr_config.xml
    Path(root, "dimr_config.xml").unlink(missing_ok=False)

    # read again, without dimr_config.xml being present
    mod2 = DFlowFMModel(root=root, mode="r")
    # check if the updated mdu was read and not newly initialized
    assert mod2.mdu.data["geometry"]["bedlevuni"] == -983


def test_write_read_model_without_geoms_crs(tmpdir):
    crs = 3857
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(root=root, mode="w", crs=crs)
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    mod1.write()

    # read again, crs now comes from geoms since it is missing in the mesh
    mod2 = DFlowFMModel(root=root, mode="r")
    assert mod1.crs.to_epsg() == crs
    assert mod2.crs.to_epsg() == crs
    # if the user provides a different crs, the correct crs still comes from the geoms
    # but only after reading the mesh/model
    mod3 = DFlowFMModel(root=root, mode="r", crs=4326)
    assert mod3.crs.to_epsg() == 4326
    mod3.read()
    assert mod3.crs.to_epsg() == crs

    # remove geoms and read again, it fails since the mesh does not have a crs
    shutil.rmtree(join(root, "geoms"))
    with pytest.raises(ValueError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "CRS was not found in the mesh or the geoms of the model" in str(e.value)

    # if the CRS cannot be found in the geoms/mesh, the user provided crs is set
    # this might be wrong like in this example.
    mod4 = DFlowFMModel(root=root, mode="r", crs=4326)
    assert mod4.crs.to_epsg() == 4326

    # copy the entire model to avoid renaming conflicts
    # add a crs to the mesh and read the model again
    root_new = join(tmpdir, "dflowfm_example_copy")
    shutil.copytree(root, root_new)
    netfile = join(root, "dflowfm/fm_net.nc")
    netfile_new = join(root_new, "dflowfm/fm_net.nc")
    uds = xu.open_dataset(netfile)
    uds.ugrid.set_crs(crs)
    uds.ugrid.to_netcdf(netfile_new)
    uds.close()
    mod5 = DFlowFMModel(root=root_new, mode="r")
    assert mod5.crs.to_epsg() == crs
    # if the user provides a different crs, the correct crs still comes from the mesh
    # but only after reading the mesh/model
    mod6 = DFlowFMModel(root=root_new, mode="r", crs=4326)
    assert mod6.crs.to_epsg() == 4326
    mod6.read()
    assert mod6.crs.to_epsg() == crs


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
    # TODO: the model crs is actually 32647, but is overwritten here
    # https://github.com/Deltares/hydromt_delft3dfm/issues/119
    root = join(EXAMPLEDIR, "dflowfm_local")
    model3 = DFlowFMModel(root=root, mode="r", crs=4326)
    model4 = DFlowFMModel(root=root, mode="r+", crs=4326)
    assert model3.crs.to_epsg() == 4326
    assert model4.crs.to_epsg() == 4326


def test_read_write_config_empty_paths(tmpdir):
    # Instantiate an empty model
    root = join(tmpdir, "dflowfm_example")
    model1 = DFlowFMModel(root=root, mode="w", crs=3857)
    model1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
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