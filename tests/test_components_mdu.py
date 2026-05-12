import pytest
from hydromt_delft3dfm import DFlowFMModel
from hydromt_delft3dfm.components.mdu import MDUComponent
import datetime as dt
from os.path import join

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


def test_get_model_time_invalid_timeformat():
    # it is required to pass a DFlowFMModel to MDUComponent, so initialize an empty one
    root = "./dflowfm_example"
    model = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mdu = MDUComponent(model=model)

    # invalid formats, should be %Y%m%d or %Y%m%d%H%M%S
    data = {
        "time.refdate": "2010-02-02",
        "time.startdatetime": "2010-02-02",
        "time.stopdatetime": "2010-02-06",
    }
    mdu.update(data)
    with pytest.raises(ValueError) as e:
        _, _ = mdu.get_model_time()
    assert "does not match format '%Y%m%d%H%M%S'" in str(e.value)


def test_mdu_write_invalid_timeformat(tmpdir):
    # it is required to pass a DFlowFMModel to MDUComponent, so initialize an empty one
    root = join(tmpdir, "dflowfm_example")
    model = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    mdu = MDUComponent(model=model)

    # invalid formats, should be %Y%m%d or %Y%m%d%H%M%S
    data = {
        "time.refdate": "2010-02-02",
        "time.startdatetime": "2010-02-02",
        "time.stopdatetime": "2010-02-06",
    }
    mdu.update(data)
    from pydantic.v1.error_wrappers import ValidationError
    with pytest.raises(ValidationError) as e:
        mdu.write()
    assert "value is not a valid integer" in str(e.value)
    assert "Invalid datetime string for startDateTime" in str(e.value)
    assert "Invalid datetime string for stopDateTime" in str(e.value)
