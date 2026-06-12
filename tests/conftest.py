from os.path import abspath, dirname, join
from hydromt_delft3dfm import DFlowFMModel
import pytest

EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


@pytest.fixture
def dflowfm_2dmodel_with_localdata(tmpdir):
    model = DFlowFMModel(
        root=str(tmpdir),
        data_libs=join(EXAMPLEDIR, "data", "data_catalog_local.yaml"),
        crs=3857,
        mode="w",
    )

    model.setup_config(
        **{
            "time.startdatetime": "20200101",
            "time.stopdatetime": "20200102",
        }
    )

    # Set a small default mesh to speed up the test, since we only want to test
    # setup_spatial_uniform_meteo and not mesh setup here.
    model.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=5000,
    )

    return model