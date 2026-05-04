import pytest
from hydromt_delft3dfm import DFlowFMModel
from hydromt_delft3dfm.components.mdu import MDUComponent
import datetime as dt

def test_get_model_time_from_refdate_tunit_start_tstop():
    # it is required to pass a DFlowFMModel to MDUComponent, so initialize an empty one
    root = "./dflowfm_example"
    model = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mdu = MDUComponent(model=model)
    data = {
        "time.refdate": "20100202",
        "time.tunit": "X",
        "time.tstart": 0,
        "time.tstop": 24,
    }
    # model.setup_config(**data)
    mdu.update(data)

    with pytest.raises(ValueError) as e:
        _ = mdu.get_model_time()
    assert "tunit='X' not supported by get_model_time()" in str(e.value)

    mdu.update(data={"time.tunit": "D"})
    _, tstop = mdu.get_model_time()
    assert tstop == dt.datetime(2010, 2, 26, 0, 0, 0)

    mdu.update(data={"time.tunit": "H"})
    _, tstop = mdu.get_model_time()
    assert tstop == dt.datetime(2010, 2, 3, 0, 0, 0)

    mdu.update(data={"time.tunit": "M"})
    _, tstop = mdu.get_model_time()
    assert tstop == dt.datetime(2010, 2, 2, 0, 24, 0)

    mdu.update(data={"time.tunit": "S"})
    _, tstop = mdu.get_model_time()
    assert tstop == dt.datetime(2010, 2, 2, 0, 0, 24)


def test_get_model_time_from_startdatetime_stopdatetime():
    # it is required to pass a DFlowFMModel to MDUComponent, so initialize an empty one
    root = "./dflowfm_example"
    model = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mdu = MDUComponent(model=model)

    # passing complete datetime string
    data = {
        "time.startdatetime": "20100202120000",
        "time.stopdatetime": "20100203120000",
    }
    # model.setup_config(**data)
    mdu.update(data)
    tstart, tstop = mdu.get_model_time()
    assert tstart == dt.datetime(2010, 2, 2, 12, 0, 0)
    assert tstop == dt.datetime(2010, 2, 3, 12, 0, 0)

    # omitting hhmmss
    data = {
        "time.startdatetime": "20100202",
        "time.stopdatetime": "20100203",
    }
    # model.setup_config(**data)
    mdu.update(data)
    tstart, tstop = mdu.get_model_time()
    assert tstart == dt.datetime(2010, 2, 2, 0, 0, 0)
    assert tstop == dt.datetime(2010, 2, 3, 0, 0, 0)
