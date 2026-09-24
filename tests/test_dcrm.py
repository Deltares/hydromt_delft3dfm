from os.path import abspath, dirname, join

from hydromt_delft3dfm import DFlowFMModel

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")
TOLERANCE = 1e-6


def test_dcrm():
    """
    Writing a model for dcrm
    """
    model = DFlowFMModel(
        root="c:\\Projects\\Hydromt-Delft3dfm\\vanSocorro\\model_dcrm",
        data_libs="c:\\Projects\\Hydromt-Delft3dfm\\vanSocorro\\data_catalog_gauges_dcrm.yaml",
        crs=3857,
        mode="w",
    )

    model.setup_config(
        **{
            "geometry.bedlevuni": -5,
            "time.startdatetime": "20100101",
            "time.stopdatetime": "20100201",
        }
    )

    model.setup_mesh2d(
        region={'geom': "c:\\Projects\\Hydromt-Delft3dfm\\vanSocorro\\tekong_reservoir_3857.geojson"},
        res=50,
    )

    model.setup_timeseries_meteo(
        meteo_timeseries_fn="rainfall",
        meteo_type="rainfall",
        fill_value=0.0,
        )
    
    model.write()

