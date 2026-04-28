#!/usr/bin/env python3
"""
Dry deposition box model helper: computes the Hg0 deposition velocity
and the corresponding first-order loss rate constant for a well-mixed box.

    k [s⁻¹] = v_d [cm/s] × 0.01 / H [m]

Import ``compute_drydep_rate`` in other scripts to get the loss rate.
"""

from pathlib import Path

import numpy as np
import xarray as xr

from .drydep_functions import DEPVEL, METERO

# ── Olson land-type lookup tables ─────────────────────────────────────────────
_data_dir = Path(__file__).parent.parent / "data"
_ds = xr.open_dataset(_data_dir / "Olson_2001_Drydep_Inputs.nc")

DRYCOEFF = _ds.DRYCOEFF.values
IOLSON = _ds.IOLSON.values
IDEP = _ds.IDEP.values
IRI = _ds.IRI.values
IRI[2] = 200.0  # match current GEOS-Chem
IRLU = _ds.IRLU.values
IRAC = _ds.IRAC.values
IRGSS = _ds.IRGSS.values
IRGSO = _ds.IRGSO.values
IRCLS = _ds.IRCLS.values
IRCLO = _ds.IRCLO.values

# ── Default met / land-cover conditions (tropical mixed-forest) ───────────────
DEFAULT_MET = dict(
    TC0=302.0,  # surface temperature (K)
    CFRAC=0.8,  # cloud fraction
    RADIAT=606.0,  # incident shortwave (W m⁻²)
    AZO=1.5,  # roughness height (m)
    USTAR=0.15,  # friction velocity (m s⁻¹)
    PRESSU=1e5,  # surface pressure (Pa)
    SUNCOS_MID=0.8,  # cos(solar zenith angle)
    ALBD=0.1,  # albedo
    BXHEIGHT=130.0,  # reference box height for METERO call (m)
    U10M=0.1,  # zonal 10-m wind (m s⁻¹)
    V10M=-0.05,  # meridional 10-m wind (m s⁻¹)
    AIRDEN=1.15,  # dry air density (kg m⁻³)
    HFLUX=16.8,  # sensible heat flux (W m⁻²)
)

_IREG = 4
_ILAND = np.zeros(73, dtype=int)
_ILAND[0] = 0
_ILAND[1] = 29
_ILAND[2] = 33
_ILAND[3] = 43
_IUSE = np.zeros(73)
_IUSE[0] = 312.0
_IUSE[1] = 200.0
_IUSE[2] = 175.0
_IUSE[3] = 313.0
_XLAI = np.zeros(73)
_XLAI[0] = 0.1
_XLAI[29] = 0.71
_XLAI[33] = 0.61
_XLAI[43] = 0.74

# Hg0 species parameters
XMW = 201e-3  # molar mass (kg mol⁻¹)
HSTAR = 0.11  # Henry's Law constant (M atm⁻¹)
F0 = 3e-5  # reactivity factor


def compute_drydep_rate(box_height_m=1.0, met=None):
    """Return the Hg0 dry deposition velocity and first-order loss rate.

    Parameters
    ----------
    box_height_m : float
        Height of the well-mixed box (m).
    met : dict, optional
        Meteorological override dict; missing keys fall back to DEFAULT_MET.

    Returns
    -------
    dvel_cms : float
        Deposition velocity (cm s⁻¹).
    k_s : float
        First-order loss rate constant (s⁻¹).
    """
    m = DEFAULT_MET if met is None else {**DEFAULT_MET, **met}

    CZ1, LSNOW, OBK, _ = METERO(
        m["BXHEIGHT"],
        m["ALBD"],
        m["TC0"],
        m["USTAR"],
        m["AIRDEN"],
        m["HFLUX"],
        m["U10M"],
        m["V10M"],
    )

    dvel_cms = DEPVEL(
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
        _IREG,
        _ILAND,
        _IUSE,
        m["TC0"],
        _XLAI,
        LSNOW,
        m["RADIAT"],
        m["CFRAC"],
        m["SUNCOS_MID"],
        m["PRESSU"],
        m["USTAR"],
        m["AZO"],
        CZ1,
        OBK,
        XMW,
        F0,
        HSTAR,
    )

    k_s = dvel_cms * 0.01 / box_height_m
    return dvel_cms, k_s


if __name__ == "__main__":
    dvel, k = compute_drydep_rate(box_height_m=1.0)
    print(f"v_d = {dvel:.4f} cm/s,  k = {k:.4e} s⁻¹,  lifetime = {1 / k:.1f} s")
