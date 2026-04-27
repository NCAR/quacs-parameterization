#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated: Jan 2022
Running offline version of dry deposition of GEOS-Chem for simple case
@author: arifeinberg
"""
from pathlib import Path

import numpy as np
import xarray as xr

from quacs.drydep.simple.drydep_functions import DEPVEL, METERO

# ── Load necessary parameters for running dry deposition code ─────────────────

data_dir = Path(__file__).parent.parent.parent / "data"
ds_ols1 = xr.open_dataset(data_dir / "Olson_2001_Drydep_Inputs.nc")

DRYCOEFF = ds_ols1.DRYCOEFF.values
IOLSON = ds_ols1.IOLSON.values
IDRYDEP = ds_ols1.IDRYDEP.values

IDEP = ds_ols1.IDEP.values
IZO = ds_ols1.IZO.values
IRI = ds_ols1.IRI.values
IRI[2] = 200.0  # match current GEOS-Chem
IRLU = ds_ols1.IRLU.values
IRAC = ds_ols1.IRAC.values
IRGSS = ds_ols1.IRGSS.values
IRGSO = ds_ols1.IRGSO.values
IRCLS = ds_ols1.IRCLS.values
IRCLO = ds_ols1.IRCLO.values
IVSMAX = ds_ols1.IVSMAX.values

# ── Inputs to the dry deposition model ────────────────────────────────────────

TC0 = 302  # temperature (K)
CFRAC = 0.8  # cloud fraction (unitless)
RADIAT = 606  # incident shortwave at ground (W/m2)
AZO = 1.5  # roughness height (m)
USTAR = 0.15  # friction velocity (m/s)
PRESSU = 1e5  # surface pressure (Pa)
SUNCOS_MID = 0.8  # cosine of solar zenith angle (unitless)
ALBD = 0.1  # albedo (unitless)
BXHEIGHT = 130  # height of lowermost grid box in GEOS-Chem (m)
U10M = 0.1  # zonal wind speed at 10 m height (m/s)
V10M = -0.05  # meridional wind speed at 10 m (m/s)
AIRDEN = 1.15  # dry air density (kg/m3)
HFLUX = 16.8  # sensible heat flux (W/m2)

IREG = 4  # number of land cover categories in this grid box

# Land cover category numbers refer to Olson land types, see:
# http://wiki.seas.harvard.edu/geos-chem/index.php/Olson_land_map
ILAND = np.zeros(73, dtype=int)
ILAND[0] = 0  # Water
ILAND[1] = 29  # Seasonal Tropical Forest
ILAND[2] = 33  # Tropical Rainforest
ILAND[3] = 43  # Savanna (Woods)

IUSE = np.zeros(73)
IUSE[0] = 312.0  # Water
IUSE[1] = 200.0  # Seasonal Tropical Forest
IUSE[2] = 175.0  # Tropical Rainforest
IUSE[3] = 313.0  # Savanna (Woods)

XLAI = np.zeros(73)
XLAI[0] = 0.1  # Water
XLAI[29] = 0.71  # Seasonal Tropical Forest
XLAI[33] = 0.61  # Tropical Rainforest
XLAI[43] = 0.74  # Savanna (Woods)

XMW = 201e-3  # Hg0 molar mass (kg/mol)
HSTAR = 0.11  # Hg0 Henry's Law Constant (M/atm)
F0 = 3e-5  # Hg0 reactivity

# ── Run dep velocity function ──────────────────────────────────────────────────

CZ1, LSNOW, OBK, W10 = METERO(BXHEIGHT, ALBD, TC0, USTAR, AIRDEN, HFLUX, U10M, V10M)

DV = DEPVEL(
    DRYCOEFF,
    IOLSON,
    IDEP,
    IRI,
    IRLU,
    IRAC,
    IRGSS,
    IRGSO,
    IRCLS,
    IRCLO,
    IREG,
    ILAND,
    IUSE,
    TC0,
    XLAI,
    LSNOW,
    RADIAT,
    CFRAC,
    SUNCOS_MID,
    PRESSU,
    USTAR,
    AZO,
    CZ1,
    OBK,
    XMW,
    F0,
    HSTAR,
)

print("Calculated dry deposition velocity (cm/s): ")
print(DV)
