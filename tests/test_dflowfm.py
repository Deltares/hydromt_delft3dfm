from os.path import abspath, basename, dirname, join
from os import makedirs, rename
from hydromt.data_catalog import DataCatalog
from hydromt.error import NoDataException
from hydromt_delft3dfm import DFlowFMModel
import numpy as np
from pathlib import Path
import pytest
import shutil
import xugrid as xu

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")
TOLERANCE = 1e-6


def _write_csv(outputdir, lines):
    meteo_fn = join(outputdir, "meteo_timeseries.csv")
    with open(meteo_fn, "w") as f:
        f.write("\n".join(lines))


def _model_update_datacatalog(model, datacat_contents):
    model_root = model.root.path
    datacat_file = join(model_root, "dummy_data_catalog.yaml")
    with open(datacat_file, "w") as f:
        f.write(datacat_contents)
    datacat = DataCatalog(datacat_file)
    model.data_catalog.update_sources(
        meteo_timeseries=datacat.get_source("meteo_timeseries"),
    )


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
    assert f"{basename(root)} does not exist" in str(e.value)

    # create the root directory
    makedirs(root, exist_ok=False)
    # point to an empty folder fails because the mdu file cannot be found.
    with pytest.raises(FileNotFoundError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "hydromt_delft3dfm requires an mdu file in read mode" in str(e.value)

    # create a dummy mdu file to move to the next error
    mdu_filename = "dflowfm/DFlowFM.mdu"
    mdu_filepath = join(root, mdu_filename)
    makedirs(dirname(mdu_filepath), exist_ok=False)
    with open(mdu_filepath, "w") as f:
        f.write("")
    # point to an empty folder with only an mdu file fails because the network file
    # cannot be found. This only happens because crs=None, so DFlowFMModel._check_crs()
    # activates the reading of the mesh during the initialisation of the DFlowFMModel.
    with pytest.raises(ValueError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "hydromt_delft3dfm cannot read a model without a mesh/network" in str(e.value)

    # when providing a crs, DFlowFMModel._check_crs() does not activate the reading the
    # mesh during initialisation of the DFlowFMModel, but then the reading fails later
    # when calling something that tries to read the mesh.
    mod3 = DFlowFMModel(root=root, mode="r", crs=4326)
    with pytest.raises(ValueError) as e:
        mod3.read()
    assert "hydromt_delft3dfm cannot read a model without a mesh/network." in str(e.value)


def test_write_read_mesh_model_different_dimr_mdu_mesh_path(tmpdir):
    """
    A model consists of at least a dimr (optional on read), a mdu and a network.
    However, the paths can be different for each model, so check if write/read
    still works when using non-default values.
    And check if the crs is preserved with write/read.
    """
    crs = 3857
    fn_dimr = "dimr.xml"
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=crs,
        mdu_filename="folder/nonstandard.mdu",
        dimr_filename=fn_dimr,
    )
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
    # check if the dimr file has the correct fm paths
    # components are only available after write
    assert str(mod1.dimr.data.component[0].workingDir) == 'folder'
    assert str(mod1.dimr.data.component[0].inputFile) == 'nonstandard.mdu'
    # check if the mdu keywords were indeed updated
    assert mod1.mdu.data["geometry"]["bedlevuni"] == geom_bedlevuni
    assert str(mod1.mdu.data["geometry"]["netfile"]) == geom_netfile
    # check if the network with the non-default path was found
    assert mod1.mesh.is_empty is False

    # read in the model to see if all changes are preserved
    mod2 = DFlowFMModel(root=root, mode="r", dimr_filename=fn_dimr)
    assert mod1.crs.to_epsg() == crs
    assert mod2.crs.to_epsg() == crs

    # check if the dimr file has the correct fm paths
    assert str(mod2.dimr.data.component[0].workingDir) == 'folder'
    assert str(mod2.dimr.data.component[0].inputFile) == 'nonstandard.mdu'
    # check if the updated mdu was read and not newly initialized
    assert mod2.mdu.data["geometry"]["bedlevuni"] == geom_bedlevuni
    assert str(mod2.mdu.data["geometry"]["netfile"]) == geom_netfile
    # check if the network with the non-default path was found
    assert mod2.mesh.is_empty is False


def test_write_read_model_without_dimr_mdu(tmpdir):
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

    # also remove the mdu file
    Path(root, "dflowfm/DFlowFM.mdu").unlink(missing_ok=False)

    # read again, this fails without dimr_config.xml and mdu
    with pytest.raises(FileNotFoundError) as e:
        _ = DFlowFMModel(root=root, mode="r")
    assert "hydromt_delft3dfm requires an mdu file in read mode" in str(e.value)


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
    assert model1.mdu.data["output"]["waqoutputdir"] == ""
    assert model1.mdu.data["trachytopes"]["trtdef"] == ""
    assert model1.mdu.data["trachytopes"]["trtl"] == ""

    # write the model to read it again
    model1.write()
    model2 = DFlowFMModel(root=root, mode="r", crs=3857)
    # Get the mdu settings
    model2.mdu.read()
    # Check whether the path is an emtpy string (was Path("") before)
    # was fixed in https://github.com/Deltares/HYDROLIB-core/issues/703
    # and https://github.com/Deltares/HYDROLIB-core/issues/1053
    assert model2.mdu.data["output"]["outputdir"] == ""
    assert model2.mdu.data["output"]["waqoutputdir"] == ""
    assert model2.mdu.data["trachytopes"]["trtdef"] == ""
    assert model2.mdu.data["trachytopes"]["trtl"] == ""


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


def test_setup_rivers_from_dem(tmpdir):
    """
    based on test_model_build[piave]
    also raises the NumbaTypeSafetyWarning to be resolved in
    https://github.com/Deltares/hydromt_delft3dfm/issues/289
    """
    root = join(tmpdir, "dflowfm_example")
    model = DFlowFMModel(
        root=root, mode="w", crs=3857, data_libs=["artifact_data"],
    )
    model.setup_rivers_from_dem(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        hydrography_fn="merit_hydro",
        river_geom_fn="hydro_rivers_lin",
    )
    # the river width values changed after a bugfix in pyflwdir, more info in
    # https://github.com/Deltares/hydromt_delft3dfm/issues/297
    rivwidth_actual = model.geoms.data["rivers"]["width"].values
    # rivwidth_expected = np.array(
    #     [60.50078133, 63.595319  , 50.        , 60.93007025, 60.93007025,
    #      ]) # pyflwdir pypi
    rivwidth_expected = np.array(
        [55.000397  , 64.28401   , 50.        , 60.93007025, 60.93007025,
         ]) # pyflwdir main
    assert np.allclose(rivwidth_actual, rivwidth_expected)


def test_setup_channels(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")

    # setup_channels
    region = {'geom': join(EXAMPLEDIR, "data", "local_data","1D_extent.geojson")}
    channels_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers.geojson")
    crosssections_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers_pointcrosssections.geojson")
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
    retentions_fn = join(EXAMPLEDIR, "data", "local_data","retention_ponds.geojson")
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
    region = {'geom': join(EXAMPLEDIR, "data", "local_data","1D_extent.geojson")}
    channels_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers.geojson")
    crosssections_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )

    # setup bridges (total of 2 bridges)
    bridges_fn = join(EXAMPLEDIR, "data", "local_data","bridges.geojson")
    model.setup_bridges(bridges_fn=bridges_fn)
    assert len(model.geoms.data['bridges']) == 2 


def test_setup_culverts(tmpdir):
    # Instantiate a dummy model
    model = DFlowFMModel(root=join(EXAMPLEDIR, "dflowfm_local"), mode="r")
    model.read()
    model.root.set(tmpdir, mode="w")

    # first add channels to obtain friction values for branches
    # see also https://github.com/Deltares/hydromt_delft3dfm/issues/168
    region = {'geom': join(EXAMPLEDIR, "data", "local_data","1D_extent.geojson")}
    channels_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers.geojson")
    crosssections_fn = join(EXAMPLEDIR, "data", "local_data","1D_rivers_pointcrosssections.geojson")
    model.setup_channels(
        region=region, channels_fn=channels_fn,
        crosssections_fn=crosssections_fn,
        crosssections_type='point'
    )

    # setup culverts (total of 1 culvert)
    culverts_fn = join(EXAMPLEDIR, "data", "local_data","culverts.geojson")
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
    raster_fn = join(EXAMPLEDIR, "data", "local_data","frictioncoefficient.tif")
    variable = 'roughness_manning'
    variables = [variable]
    model.inifield.add_raster_data_from_rasterdataset(raster_fn, variables)

    roughness_values = np.unique(model.inifield.data[variable]).tolist()
    expected_values = [-999.0, 0.025, 0.044, 0.050, 0.055]
    assert np.allclose(roughness_values, expected_values, atol=TOLERANCE)


def test_setup_spatial_forcing(tmpdir):
    root = join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        data_libs=["artifact_data"],
        crs=3857,
    )

    # change the start/stop times to be in the period of the artifact_data
    mod1.setup_config(**{
        "time.startdatetime": "20100202",
        "time.stopdatetime": "20100203",
    })

    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    # Possible variable names can be found in the translate_utils module
    #  hydromt_delft3dfm.utils.translate_utils.
    # Beware of unit conversions: https://github.com/Deltares/hydromt_delft3dfm/issues/304
    mod1.setup_spatial_forcing(
        meteo_fn="era5_hourly",  # source for meteo data
        variables=[
            "precip", "press_msl",
            # more variables available in deltares_data
            # "temp_dew", "wind10_u", "wind10_v",
            # even more in earthdatahub_data
            # "u10n", "v10n", "chnk",
        ],
    )

    # write calls the validators and writes all the delft3dfm files
    mod1.write()

    # check if the created model with netcdf forcing can also be read properly
    mod2 = DFlowFMModel(
        root=root,
        mode="r",
    )
    mod2.read()
    expected_keys = set(['rainfall', 'airpressure'])
    assert set(mod2.forcing.data.keys()) == expected_keys


def test_setup_constant_meteo(dflowfm_2dmodel_with_localdata):
    dflowfm_2dmodel_with_localdata.setup_constant_meteo(
        meteo_type="rainfall",
        constant_value=5.0,
    )
    assert "meteo_rainfall" in dflowfm_2dmodel_with_localdata.forcing.data
    mdu_rainfaill = dflowfm_2dmodel_with_localdata.mdu.get_value(
        'external_forcing.rainfall'
    )
    assert mdu_rainfaill == 1
    # to at least call the writer in one of the tests
    dflowfm_2dmodel_with_localdata.forcing.write()


def test_setup_timeseries_rainfall_rate_from_datacatalog(dflowfm_2dmodel_with_localdata):
    dflowfm_2dmodel_with_localdata.setup_timeseries_meteo(
        meteo_type="rainfall_rate",
        meteo_timeseries_fn="meteo_timeseries_T2",
    )

    assert "meteo_rainfall_rate" in dflowfm_2dmodel_with_localdata.forcing.data

    # to at least call the writer in one of the tests
    dflowfm_2dmodel_with_localdata.forcing.write()


def test_setup_timeseries_rainfall_timeseries_fills_missing_values(
        dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path
    # write meteo timeseries with missing values (nan)
    _write_csv(
        model_root,
        [
            "time,rainfall",
            "2020-01-01 00:00,2.0",
            "2020-01-02 00:00,NaN",
        ],
    )

    dflowfm_2dmodel_empty.setup_timeseries_meteo(
        meteo_type="rainfall",
        meteo_timeseries_fn="meteo_timeseries",
        fill_value=0.0,
    )

    da = dflowfm_2dmodel_empty.forcing.data["meteo_rainfall"]

    assert not np.isnan(da.values).any()
    assert np.isclose(da.values[0, -1], 0.0)


def test_setup_constant_meteo_rejects_unknown_type(dflowfm_2dmodel_with_localdata):
    with pytest.raises(ValueError, match="Unsupported meteo_type"):
        dflowfm_2dmodel_with_localdata.setup_constant_meteo(
            meteo_type="evapotranspiration",
            constant_value=1.0,
       )


def test_setup_timeseries_meteo_no_csv(
     dflowfm_2dmodel_empty,
):
    with pytest.raises(NoDataException, match="Resolver 'convention' found no files"):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_rejects_single_timestep(
     dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path
    # create a meteo timeseries with a single timeseries to trigger the error
    _write_csv(
        model_root,
        [
            "time,rainfall",
            "2020-01-01 00:00,2.0",
        ],
    )

    with pytest.raises(ValueError, match="must contain at least two timesteps"):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_too_short_timeseries(
     caplog,
    dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path
    # create a meteo timeseries with a single timeseries to trigger the error
    _write_csv(
        model_root,
        [
            "time,rainfall",
            "2019-12-31 00:00,2.0",
            "2019-12-31 12:00,2.0",
            "2020-01-01 00:00,2.0",
            "2020-01-01 12:00,2.0",
        ],
    )

    # with pytest.raises(ValueError, match="must contain at least two timesteps"):
    dflowfm_2dmodel_empty.setup_timeseries_meteo(
        meteo_type="rainfall",
        meteo_timeseries_fn="meteo_timeseries",
    )

    # hydromt-core warning that timeseries is too long and will be clipped
    assert "Requested time range" in caplog.text
    assert "partially overlaps with available range" in caplog.text
    assert "Clamping to (2020-01-01 00:00:00, 2020-01-01 12:00:00)" in caplog.text
    # hydromt_delft3dfm warning that timeseries will be padded with fill_value
    assert "Time in meteo_timeseries_fn is shorter than the model" in caplog.text
    assert "Missing values will be filled using 0.0" in caplog.text

    # assert resulting timeseries
    ts = dflowfm_2dmodel_empty.forcing.data["meteo_rainfall"].to_numpy()
    assert np.allclose(ts, [[2., 2., 0.]])


def test_setup_timeseries_meteo_rejects_non_equidistant_timeseries(
     dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path
    # create a non-equidistant meteo timeseries to trigger the error
    # meteo_timeseries.csv is predefined in the data_catalog.yaml in the
    # dflowfm_2dmodel_empty fixture
    _write_csv(
        model_root,
        [
            "time,rainfall",
            "2020-01-01 00:00,2.0",
            "2020-01-01 00:20,2.0",
            "2020-01-02 00:00,2.0",
        ],
    )

    with pytest.raises(ValueError, match="Non-equidistant time series"):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_rejects_unknown_freq(
     dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path
    # create a ts with ms frequency to trigger the error
    _write_csv(
        model_root,
        [
            "time,rainfall",
            "2020-01-01 00:00:00.000,2.0",
            "2020-01-01 00:00:00.200,2.0",
            "2020-01-02 00:00:01.400,2.0",
        ],
    )

    with pytest.raises(ValueError, match="Unsupported time frequency 'ms'"):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_rejects_no_column_named_time(
    dflowfm_2dmodel_empty,
):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path

    # update the meteo_timeseries datasource to have index_col=time (instead of 0),
    #  but the csv deliberately has a column called `date`, triggering the error.
    # TODO: even if we would provide the correct column name (date) there will still be
    #  an error until https://github.com/Deltares/hydromt/issues/1502 is fixed.
    #  The only way to get it working at the moment is with index_col=0.
    datacat_contents = """
        meteo_timeseries:
          data_type: DataFrame
          uri: meteo_timeseries.csv
          driver:
            name: pandas
            options:
              index_col: time
              parse_dates: true
          metadata:
            unit: mm day-1
        """
    _model_update_datacatalog(dflowfm_2dmodel_empty, datacat_contents)

    _write_csv(
        model_root,
        [
            "date,rainfall",
            "2020-01-01 00:00,2.0",
            "2020-01-02 00:00,2.0",
        ],
    )

    # The error is "'time' not in list" in python<=3.13, but this has changed to
    #  "list.index(x): x not in list" in python 3.14. This can be reproduced with:
    #  `["a", "b", "c"].index("d")`. Therefore only match the end of the error message.
    with pytest.raises(ValueError, match=" not in list"):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_rejects_no_time_index_from_datacatalog(dflowfm_2dmodel_empty):
    # the model_root is a tmpdir named to the test that calls the fixture
    model_root = dflowfm_2dmodel_empty.root.path

    # create dummy catalog with incomplete driver (commented)
    # this test is purely to trigger the error
    datacat_contents = """
meteo_timeseries:
  data_type: DataFrame
  uri: meteo_timeseries.csv
  driver:
    name: pandas
    #options:
    #  index_col: 0
    #  parse_dates: true
  metadata:
    unit: mm day-1
"""
    _model_update_datacatalog(dflowfm_2dmodel_empty, datacat_contents)

    _write_csv(
        model_root,
        [
            "time,rainfall_rate",
            "2020-01-01 00:00,2.0",
            "2020-01-02 00:00,2.0",
        ],
    )

    err_msg = "meteo_timeseries_fn must provide a datetime index"
    with pytest.raises(ValueError, match=err_msg):
        dflowfm_2dmodel_empty.setup_timeseries_meteo(
            meteo_type="rainfall_rate",
            meteo_timeseries_fn="meteo_timeseries",
        )


def test_setup_timeseries_meteo_rejects_no_matching_variable(
    dflowfm_2dmodel_with_localdata,
):
    with pytest.raises(ValueError, match="columns expected but not found"):
        dflowfm_2dmodel_with_localdata.setup_timeseries_meteo(
            meteo_type="rainfall",
            meteo_timeseries_fn="meteo_timeseries_T2",
        )


def test_setup_rainfall_from_constant_deprecated(
    dflowfm_2dmodel_with_localdata,
):
    err_msg = "setup_rainfall_from_constant is deprecated"
    with pytest.raises(AttributeError, match=err_msg):
        dflowfm_2dmodel_with_localdata.setup_rainfall_from_constant(constant_value=5.0)


def test_setup_rainfall_from_uniform_timeseries_deprecated(
    dflowfm_2dmodel_with_localdata,
):
    err_msg = "setup_rainfall_from_uniform_timeseries is deprecated"
    with pytest.raises(AttributeError, match=err_msg):
        dflowfm_2dmodel_with_localdata.setup_rainfall_from_uniform_timeseries(
            meteo_timeseries_fn="meteo_timeseries_T2",
            fill_value=0.0,
            is_rate=True,
        )
