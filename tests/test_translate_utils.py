import pytest
from hydromt_delft3dfm.utils.translate_utils import varname_to_dflowfm_quantity

def test_varname_to_dflowfm_quantity_hydromt():
    varname = "precip"
    quantity = varname_to_dflowfm_quantity(varname)
    assert quantity == "rainfall"


def test_varname_to_dflowfm_quantity_era5():
    varname = "u10"
    quantity = varname_to_dflowfm_quantity(varname)
    assert quantity == "windx"


def test_varname_to_dflowfm_quantity_quantity():
    varname = "airpressure"
    quantity = varname_to_dflowfm_quantity(varname)
    assert quantity == "airpressure"


def test_varname_to_dflowfm_quantity_notpresent():
    varname = "u10_typo"
    with pytest.raises(KeyError) as e:
        _ = varname_to_dflowfm_quantity(varname)
    assert "varname u10_typo not in keys or values of translation" in str(e.value)
