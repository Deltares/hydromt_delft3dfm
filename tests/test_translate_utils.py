import pytest
from hydromt_delft3dfm.utils.translate_utils import (
    meteo_unit_from_type,
    varname_to_dflowfm_quantity,
)
from hydromt_delft3dfm.utils.translate_utils import (
    DICT_VARNAME_TO_DFLOWFM,
    METEO_UNITS,
)

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


def test_meteo_unit_from_type():
    unit = meteo_unit_from_type("airpressure")
    assert unit == "N/m2"


def test_meteo_unit_from_type_unsupported_quantity():
    with pytest.raises(ValueError) as e:
        _ = meteo_unit_from_type("aaa")
    assert "Unsupported meteo_type 'aaa'. Supported values are" in str(e.value)


def test_consistency_varnames_units():
    exceptions = ["sea_surface_temperature"]
    for quantity in DICT_VARNAME_TO_DFLOWFM.values():
        if quantity in METEO_UNITS.keys():
            continue
        if quantity in exceptions:
            continue
        raise KeyError(f"quantity '{quantity}' not in METEO_UNITS keys.")
