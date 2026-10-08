"""Wesely (1989) gas dry deposition, translated from CAM ``mo_drydep.F90``.

>>> import numpy as np
>>> from quacs.drydep.wesely import wesely_gas
>>> frac = np.zeros(11); frac[3] = 1.0   # deciduous forest
>>> out = wesely_gas("O3", sfc_temp=298., air_temp=296., pressure_sfc=101325.,
...                  pressure_10m=100000., wind_speed=5., spec_hum=0.01,
...                  solar_flux=400., frac_landuse=frac, month=7, lat=40.)
>>> out["Vd"]  # cm/s
"""

from .wesely import (
    LAND_NAMES,
    N_LANDUSE,
    N_SEASONS,
    SPECIES_TABLE,
    GasSpecies,
    calc_heff,
    calc_ra_rb,
    calc_rc,
    detect_dew,
    get_season_index,
    run_box_model,
    wesely_gas,
)

__all__ = [
    "N_LANDUSE",
    "N_SEASONS",
    "LAND_NAMES",
    "GasSpecies",
    "SPECIES_TABLE",
    "calc_heff",
    "get_season_index",
    "calc_ra_rb",
    "calc_rc",
    "detect_dew",
    "wesely_gas",
    "run_box_model",
]
