"""Utilities translate variable/quantity names to Delft3D FM conventions."""
import logging

logger = logging.getLogger(f"hydromt.{__name__}")

# dict to translate hydromt to dflowfm variable name conventions
DICT_VARNAME_TO_DFLOWFM = {
    # Hydromt to dflowfm (from hydromt deltares_data datacatalog
    #  era5_hourly.data_adapter.rename or https://deltares.github.io/hydromt/stable/
    #  user_guide/data_catalog/data_conventions.html#meteorology)
    "wind10_u": "windx",  # ERA5 u10
    "wind10_v": "windy",  # ERA5 v10
    "press_msl": "airpressure",  # ERA5 msl
    "temp": "airtemperature",  # ERA5 t2m
    "temp_dew": "dewpoint",  # ERA5 d2m
    "precip": "rainfall",  # ERA5 tp
    # "kin": "", # ERA5 ssrd
    # "kout": "", # ERA5 tisr
    # ERA5 to dflowfm (from dfm_tools.modelbuilder.preprocess_merge_meteofiles_era5 and
    #  long names from dfm_tools.download.download_ERA5).
    # These can be retrieved from the ERA5 data if they were not renamed in the
    #  data_catalog for instance u10n is not in the hydromt naming conventions but
    #  can in this way be retrieved from the pre-defined earthdatahub_data datacatalog.
    "u10": "windx",  # ERA5 long: 10m_u_component_of_wind
    "u10n": "windx",  # ERA5 long: 10m_u_component_of_neutral_wind
    "v10": "windy",  # ERA5 long: 10m_v_component_of_wind
    "v10n": "windy",  # ERA5 long: 10m_v_component_of_neutral_wind
    "msl": "airpressure",  # ERA5 long: mean_sea_level_pressure
    "t2m": "airtemperature",  # ERA5 long: 2m_temperature
    # TODO: when adding dewpointtemperature/netsolarradiation/others, also change
    #  mdu.physics.temperature = 5
    "d2m": "dewpoint",  # ERA5 long: 2m_dewpoint_temperature
    "ssr": "netsolarradiation",  # ERA5 long: surface_net_solar_radiation
    "tcc": "cloudiness",  # ERA5 long: total_cloud_cover
    "tp": "rainfall",
    "sst": "sea_surface_temperature",  # ERA5 long: sea_surface_temperature
    "strd": "longwaveradiation",  # ERA5 long: surface_thermal_radiation_downwards
    "chnk": "charnock",  # ERA5 long: charnock
    # TODO: mer/avg_ie is negative precipitation, also requires operand="+"
    # https://github.com/Deltares/hydromt_delft3dfm/issues/301
    "mer": "rainfall_rate",  # ERA5 long: mean_evaporation_rate
    "mtpr": "rainfall_rate",  # ERA5 long: mean_total_precipitation_rate
    # mer and mtpr are now called avg_ie and avg_tprate
    "avg_ie": "rainfall_rate",  # ERA5 long: mean_evaporation_rate
    "avg_tprate": "rainfall_rate",  # ERA5 long: mean_total_precipitation_rate
    # TODO: rhoao/mer/mtpr/avg_ie/avg_tprate are not available in earthdatahub_data
    #  probably not both mer/avg_ie and mtpr/avg_tprate will be included anyway
    #  https://github.com/Deltares/hydromt_delft3dfm/issues/294
    "rhoao": "airdensity",  # ERA5 long: air_density_over_the_oceans
    # "slhf":"",  # ERA5 long: surface_latent_heat_flux
    # "sshf":"",  # ERA5 long: surface_sensible_heat_flux
    # "str":"",  # ERA5 long: surface_net_thermal_radiation
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
