"""Utilities translate variable/quantity names to Delft3D FM conventions."""
import logging

logger = logging.getLogger(f"hydromt.{__name__}")

# dict to translate hydromt to dflowfm variable name conventions
DICT_VARNAME_TO_DFLOWFM = {
    # Hydromt to dflowfm (from hydromt deltares_data datacatalog
    #  era5_hourly.data_adapter.rename or https://deltares.github.io/hydromt/stable/
    #  user_guide/data_catalog/data_conventions.html#meteorology)
    "precip": "rainfall",  # ERA5 tp
    "wind10_u": "windx",  # ERA5 u10
    "wind10_v": "windy",  # ERA5 v10
    "press_msl": "airpressure",  # ERA5 msl
    "temp": "airtemperature",  # ERA5 t2m
    "temp_dew": "dewpoint",  # ERA5 d2m
    "kin": "solarradiation",  # ERA5 ssrd
    # ERA5 to dflowfm (from dfm_tools.modelbuilder.preprocess_merge_meteofiles_era5 and
    #  long names from dfm_tools.download.download_ERA5).
    # These can be retrieved from the ERA5 data if they were not renamed in the
    #  data_catalog for instance u10n is not in the hydromt naming conventions but
    #  can in this way be retrieved from the pre-defined earthdatahub_data datacatalog.
    "u10": "windx",  # CDS varname: 10m_u_component_of_wind
    "u10n": "windx",  # CDS varname: 10m_u_component_of_neutral_wind
    "v10": "windy",  # CDS varname: 10m_v_component_of_wind
    "v10n": "windy",  # CDS varname: 10m_v_component_of_neutral_wind
    "msl": "airpressure",  # CDS varname: mean_sea_level_pressure
    "chnk": "charnock",  # CDS varname: charnock
    "sst": "sea_surface_temperature",  # CDS varname: sea_surface_temperature
    "t2m": "airtemperature",  # CDS varname: 2m_temperature
    # TODO: when adding dewpointtemperature/netsolarradiation/others,
    #  also update `mdu.physics.temperature = 5`
    "d2m": "dewpoint",  # CDS varname: 2m_dewpoint_temperature
    "tcc": "cloudiness",  # CDS varname: total_cloud_cover
    "ssrd": "solarradiation",  # CDS varname: surface_solar_radiation_downwards
    "ssr": "netsolarradiation",  # CDS varname: surface_net_solar_radiation
    "strd": "longwaveradiation",  # CDS varname: surface_thermal_radiation_downwards
    "tp": "rainfall",  # CDS varname: total_precipitation
    # TODO: mer/avg_ie is negative precipitation, this requires operand="+"
    # https://github.com/Deltares/hydromt_delft3dfm/issues/301
    "mer": "rainfall_rate",  # CDS varname: mean_evaporation_rate
    "mtpr": "rainfall_rate",  # CDS varname: mean_total_precipitation_rate
    # mer and mtpr are now called avg_ie and avg_tprate
    "avg_ie": "rainfall_rate",  # CDS varname: mean_evaporation_rate
    "avg_tprate": "rainfall_rate",  # CDS varname: mean_total_precipitation_rate
    # TODO: rhoao/mer/mtpr/avg_ie/avg_tprate are not available in earthdatahub_data
    #  probably not both mer/avg_ie and mtpr/avg_tprate will be included anyway
    #  https://github.com/Deltares/hydromt_delft3dfm/issues/294
    "rhoao": "airdensity",  # CDS varname: air_density_over_the_oceans
}


def varname_to_dflowfm_quantity(varname: str):
    """Get quantity names from the translation dict.

    If varname is in dict.values(), quantity is equal to varname,
    if varname is in dict.keys(), quantity is the corresponding value and
    raise KeyError if varname is not in dict keys or values.

    For instance the pre-defined earthdatahub_data catalog renames some era5 variables
    to the hydromt conventions, but the other era5 variables can also be retrieved
    since the zarr archive contains all era5 variables and not only the ones renamed by
    the data catalog.
    """
    if varname in DICT_VARNAME_TO_DFLOWFM.values():
        quantity = varname
        logger.info(f"varname {varname} is already a dflowfm quantity.")
    elif varname in DICT_VARNAME_TO_DFLOWFM.keys():
        quantity = DICT_VARNAME_TO_DFLOWFM[varname]
        logger.info(f"varname {varname} translated to dflowfm quantity {quantity}.")
    else:
        raise KeyError(
            f"varname {varname} not in keys or values of translation dictionary."
        )
    return quantity
