#!/usr/bin/env python3
"""
Pythonic dry deposition physics functions.

Equivalent to the functions in ``drydep_functions.py`` but written in
idiomatic Python/NumPy:

- Snake_case names throughout
- No hard-coded array index look-up tables; resistance values are carried
  directly on ``LandCoverPatch`` objects
- ``np.inf`` instead of sentinel 9999 values (``1 / np.inf == 0`` in NumPy,
  which eliminates all downstream guard branches)
- Per-patch loop in ``deposition_velocity`` replaced by vectorised NumPy
  array operations
- Aerodynamic and quasi-laminar resistances (species/met only, patch-agnostic)
  hoisted out of the per-patch computation
- ``RA`` naming collision in the original fixed: shortwave radiation → ``radiat``
  (never reassigned); aerodynamic resistance → ``r_aero``

Public API
----------
monin_obukhov_length(ts, ustar, airden, hflux) -> float
calc_met_vars(box_height, albedo, ts, ustar, airden, hflux, u10m, v10m) -> MetVars
molecular_diffusivity(temp_k, pressure, mol_weight) -> float
deposition_velocity(patches, *, ts, radiat, cloud_frac, suncos_mid,
                    pressure, ustar, roughness_length, cz1, obk,
                    mol_weight, f0, hstar, drycoeff) -> float  [cm s⁻¹]

Internal
--------
_light_correction(drycoeff, lai, suncos, cloud_frac) -> float

References
----------
Wesely (1989) Atmos. Environ.
Wang et al. (1998) JGR
Seinfeld (1986) ch. 8
Levine (1988) ch. 15
"""

from itertools import combinations_with_replacement
from typing import NamedTuple

import numpy as np

# ── Physical constants ────────────────────────────────────────────────────────

_VON_KARMAN = 0.4  # Von Kármán constant
_CP_AIR = 1000.0  # specific heat of air at const. pressure (J kg⁻¹ K⁻¹)
_G0 = 9.80665  # gravitational acceleration (m s⁻²)
_MW_AIR = 28.8e-3  # molar mass of moist air (kg mol⁻¹)
_R = 8.3144598  # molar gas constant (J K⁻¹ mol⁻¹)
_AVO = 6.022140857e23  # Avogadro's number (mol⁻¹)
_MW_WATER = 18.016e-3  # molar mass of water (kg mol⁻¹)
_R_AIR = 1.2e-10  # radius of an air molecule (m)
_R_X = 1.5e-10  # radius of a generic trace-gas molecule (m)
_D_AIR = 0.2e-4  # thermal diffusivity of air (m² s⁻¹)


# ── MetVars namedtuple ────────────────────────────────────────────────────────


class MetVars(NamedTuple):
    """Derived meteorological quantities for the deposition scheme.

    Attributes
    ----------
    cz1 : float
        Mid-point height of the first model level (m).
    is_snow : bool
        True when the surface is covered by snow or ice (albedo > 0.4).
    obk : float
        Monin-Obukhov length (m).
    wind10 : float
        10-m wind speed magnitude (m s⁻¹).
    """

    cz1: float
    is_snow: bool
    obk: float
    wind10: float


# ── monin_obukhov_length ──────────────────────────────────────────────────────


def monin_obukhov_length(ts: float, ustar: float, airden: float, hflux: float) -> float:
    """Monin-Obukhov length (m).

    Parameters
    ----------
    ts : float
        Surface temperature (K).
    ustar : float
        Friction velocity (m s⁻¹).
    airden : float
        Air density (kg m⁻³).
    hflux : float
        Sensible heat flux (W m⁻²).  Zero is handled safely.

    Returns
    -------
    float
        Monin-Obukhov length in metres.
    """
    _hflux = hflux if hflux != 0.0 else 1e-20
    numerator = -airden * _CP_AIR * ts * ustar**3
    denominator = _VON_KARMAN * _G0 * _hflux
    return numerator / denominator


# ── calc_met_vars ─────────────────────────────────────────────────────────────


def calc_met_vars(
    box_height: float,
    albedo: float,
    ts: float,
    ustar: float,
    airden: float,
    hflux: float,
    u10m: float,
    v10m: float,
) -> MetVars:
    """Derive meteorological variables needed by the dry deposition scheme.

    Parameters
    ----------
    box_height : float
        Height of the lowest model layer (m).
    albedo : float
        Surface albedo (dimensionless, 0–1).
    ts : float
        Surface temperature (K).
    ustar : float
        Friction velocity (m s⁻¹).
    airden : float
        Air density (kg m⁻³).
    hflux : float
        Sensible heat flux (W m⁻²).
    u10m : float
        Zonal 10-m wind component (m s⁻¹).
    v10m : float
        Meridional 10-m wind component (m s⁻¹).

    Returns
    -------
    MetVars
        Named tuple with fields ``cz1``, ``is_snow``, ``obk``, ``wind10``.
    """
    cz1 = box_height / 2.0
    is_snow = bool(albedo > 0.4)
    obk = monin_obukhov_length(ts, ustar, airden, hflux)
    wind10 = float(np.sqrt(u10m**2 + v10m**2))
    return MetVars(cz1=cz1, is_snow=is_snow, obk=obk, wind10=wind10)


# ── molecular_diffusivity ─────────────────────────────────────────────────────


def molecular_diffusivity(temp_k: float, pressure: float, mol_weight: float) -> float:
    """Molecular diffusivity of a trace gas in air (m² s⁻¹).

    Follows equations 8.5 and 8.9 of Seinfeld (1986) and 15.47 of
    Levine (1988).

    Parameters
    ----------
    temp_k : float
        Temperature (K).
    pressure : float
        Pressure (Pa).
    mol_weight : float
        Molar mass of the trace gas (kg mol⁻¹).

    Returns
    -------
    float
        Diffusivity in m² s⁻¹.
    """
    air_number_density = pressure * _AVO / (_R * temp_k)  # molec m⁻³
    collision_diameter = _R_X + _R_AIR
    mass_ratio = mol_weight / _MW_AIR  # dimensionless
    mean_free_path = 1.0 / (
        np.pi * np.sqrt(1.0 + mass_ratio) * air_number_density * collision_diameter**2
    )
    mean_speed = np.sqrt(8.0 * _R * temp_k / (np.pi * mol_weight))
    return (3.0 * np.pi / 32.0) * (1.0 + mass_ratio) * mean_free_path * mean_speed


# ── _light_correction (private) ───────────────────────────────────────────────


def _light_correction(drycoeff: np.ndarray, lai: float, suncos: float, cloud_frac: float) -> float:
    """Stomatal light-correction factor (Wang et al. 1998).

    Computes a degree-3 polynomial over normalised {1, LAI, cos(SZA), f_cloud}
    using all 20 combinations-with-replacement of the four terms.

    Parameters
    ----------
    drycoeff : ndarray, shape (20,)
        Baldocchi/Wesely polynomial coefficients.
    lai : float
        Leaf area index (m² m⁻²).
    suncos : float
        Cosine of the solar zenith angle.
    cloud_frac : float
        Column cloud fraction (0–1).

    Returns
    -------
    float
        Correction factor (≥ 0.1).
    """
    # ── Normalise inputs (inline SUNPARAM) ────────────────────────────────
    _ND = np.array([55.0, 20.0, 11.0])  # scaling factors per variable
    _X0 = np.array([11.0, 1.0, 1.0])  # maxima per variable

    raw = np.array([lai, suncos, cloud_frac], dtype=float)
    raw = np.minimum(raw, _X0)
    xlow = np.where(np.arange(3) != 2, _X0 / _ND, 0.0)
    raw = np.maximum(raw, xlow)
    normed = raw / _X0

    # Four terms: constant + three normalised variables
    terms = np.concatenate([[1.0], normed])  # shape (4,)

    # Build the 20 cubic polynomial terms C(4+3-1,3) = C(6,3) = 20
    indices = list(combinations_with_replacement(range(4), 3))  # 20 triples
    idx = np.array(indices)  # (20, 3)
    realterms = np.prod(terms[idx], axis=1)  # (20,)

    result = float(np.dot(drycoeff, realterms))
    return max(result, 0.1)


# ── deposition_velocity ───────────────────────────────────────────────────────


def deposition_velocity(
    patches,
    *,
    ts: float,
    radiat: float,
    cloud_frac: float,
    suncos_mid: float,
    pressure: float,
    ustar: float,
    roughness_length: float,
    cz1: float,
    obk: float,
    mol_weight: float,
    f0: float,
    hstar: float,
    drycoeff: np.ndarray,
    is_snow: bool = False,
) -> float:
    """Dry deposition velocity (cm s⁻¹) for one grid cell.

    Implements the Wesely (1989) resistance network.  Resistance values are
    read directly from ``LandCoverPatch`` objects; there are no index look-up
    tables.

    Parameters
    ----------
    patches : list of LandCoverPatch
        Land cover types present in the grid cell.  ``fraction`` values should
        sum to 1 (or close to it).
    ts : float
        Surface temperature (K).
    radiat : float
        Incident shortwave radiation (W m⁻²).
    cloud_frac : float
        Column cloud fraction (0–1).
    suncos_mid : float
        Cosine of the solar zenith angle at the mid-point of the timestep.
    pressure : float
        Surface pressure (Pa).
    ustar : float
        Friction velocity (m s⁻¹).
    roughness_length : float
        Aerodynamic roughness length z₀ (m).
    cz1 : float
        Mid-point height of the first model level (m).
    obk : float
        Monin-Obukhov length (m).
    mol_weight : float
        Molar mass of the species (kg mol⁻¹).
    f0 : float
        Biological reactivity factor (0–1).
    hstar : float
        Effective Henry's law constant at pH 7 (M atm⁻¹).
    drycoeff : ndarray, shape (20,)
        Baldocchi/Wesely polynomial coefficients for stomatal correction.
    is_snow : bool, optional
        If True the surface is treated as snow/ice (stomata closed).

    Returns
    -------
    float
        Deposition velocity in cm s⁻¹.
    """
    n = len(patches)
    temp_c = ts - 273.15

    # ── Kinematic viscosity of air (m² s⁻¹) ──────────────────────────────
    kin_visc = 0.151e-4 * (ts / 273.15) ** 1.77

    # ── Low-temperature additional resistance term ─────────────────────────
    r_temp = 1000.0 * np.exp(-temp_c - 4.0)

    # ── Extract per-patch resistance values ───────────────────────────────
    # Use np.inf directly where sentinel ≥ 9999 (1/inf = 0 in NumPy)
    fractions = np.array([p.fraction for p in patches])
    lai = np.array([p.lai for p in patches])

    # IUSE in the original is fraction * 1000; XLAI_scale = XLAI[IOLSON] / IUSE * 1000
    # which simplifies to lai / fraction (LAI normalised to unit area).
    # Guard against zero fraction (shouldn't occur in practice).
    lai_per_frac = np.where(fractions > 0.0, lai / fractions, 0.0)

    # Raw resistance values from patches; convert ≥ 9999 → inf
    def _to_inf(arr):
        return np.where(arr >= 9999.0, np.inf, arr)

    iri_raw = _to_inf(np.array([p.iri for p in patches]))
    irlu_raw = _to_inf(np.array([p.irlu for p in patches]))
    irac_raw = _to_inf(np.array([p.irac for p in patches]))
    irgss_raw = _to_inf(np.array([p.irgss for p in patches]))
    irgso_raw = _to_inf(np.array([p.irgso for p in patches]))
    ircls_raw = _to_inf(np.array([p.ircls for p in patches]))
    irclo_raw = _to_inf(np.array([p.irclo for p in patches]))

    # ── Snow/ice: override all resistances with those from patch 0 ───────
    # In the original DEPVEL, LSNOW forces II=1 for every land type, which
    # means ALL resistance arrays are looked up at index 0 (IRI[0], IRLU[0],
    # etc.).  In our patch-direct design, patch 0 carries those values.
    if is_snow:
        snow_iri = _to_inf(np.full(n, patches[0].iri))
        snow_irlu = _to_inf(np.full(n, patches[0].irlu))
        snow_irac = _to_inf(np.full(n, patches[0].irac))
        snow_irgss = _to_inf(np.full(n, patches[0].irgss))
        snow_irgso = _to_inf(np.full(n, patches[0].irgso))
        snow_ircls = _to_inf(np.full(n, patches[0].ircls))
        snow_irclo = _to_inf(np.full(n, patches[0].irclo))
        iri_raw, irlu_raw, irac_raw = snow_iri, snow_irlu, snow_irac
        irgss_raw, irgso_raw = snow_irgss, snow_irgso
        ircls_raw, irclo_raw = snow_ircls, snow_irclo

    # ── Stomatal resistance (RI) ──────────────────────────────────────────
    # Temperature correction factor
    if 0.0 < temp_c < 40.0:
        gfact = 400.0 / temp_c / (40.0 - temp_c)
    else:
        gfact = 100.0

    # Light correction factor (per patch — LAI varies per patch)
    gfaci = np.full(n, 100.0)
    for i in range(n):
        if radiat > 0.0 and lai_per_frac[i] > 0.0:
            gfaci[i] = 1.0 / _light_correction(drycoeff, lai_per_frac[i], suncos_mid, cloud_frac)

    # Effective stomatal resistance (inf where stomata are closed)
    ri = iri_raw * gfact * gfaci  # inf * finite = inf → correct

    # ── Cuticular resistance (RLU) ────────────────────────────────────────
    # Per-leaf irlu divided by LAI gives bulk canopy value.
    # np.errstate: division by zero when lai_per_frac=0 is intentional —
    # the np.where condition selects 1e6 for those elements.
    with np.errstate(divide="ignore", invalid="ignore"):
        rlu_base = np.where(
            (irlu_raw >= np.inf) | (lai_per_frac <= 0.0),
            1e6,
            irlu_raw / lai_per_frac,
        )
    # Temperature correction: cap at 2× base value
    rlu = np.where(
        rlu_base < 1e6,
        np.minimum(rlu_base + r_temp, 2.0 * rlu_base),
        rlu_base,
    )

    # ── Canopy aerodynamic resistance (RAC) ───────────────────────────────
    rac = np.where(irac_raw < np.inf, np.maximum(irac_raw, 1.0), np.inf)

    # ── Soil/ground resistances ───────────────────────────────────────────
    rgss_base = np.maximum(irgss_raw, 1.0)
    rgss = np.where(
        rgss_base < np.inf,
        np.minimum(rgss_base + r_temp, 2.0 * rgss_base),
        np.inf,
    )

    rgso_base = np.maximum(irgso_raw, 1.0)
    rgso = np.where(
        rgso_base < np.inf,
        np.minimum(rgso_base + r_temp, 2.0 * rgso_base),
        np.inf,
    )

    # ── Lower-canopy resistance ───────────────────────────────────────────
    rcls = np.where(
        ircls_raw < np.inf,
        np.minimum(ircls_raw + r_temp, 2.0 * ircls_raw),
        np.inf,
    )
    rclo = np.where(
        irclo_raw < np.inf,
        np.minimum(irclo_raw + r_temp, 2.0 * irclo_raw),
        np.inf,
    )

    # ── Lower-canopy aerodynamic resistance (RDC) ─────────────────────────
    # Depends on radiation; constant across patches for given met.
    rdc = 100.0 * (1.0 + 1000.0 / (radiat + 10.0))

    # ── Species-dependent surface resistance components ───────────────────
    diff_h2o = molecular_diffusivity(ts, pressure, _MW_WATER)
    diff_gas = molecular_diffusivity(ts, pressure, mol_weight)

    # Effective stomatal resistance scaled for species diffusivity + solubility
    rixx = ri * (diff_h2o / diff_gas) + 1.0 / (hstar / 3000.0 + 100.0 * f0)

    # Cuticular: original checks RLU[LDT] < 9999 (the LAI-scaled value)
    # before computing RLUXX; when rlu >= 9999, RLUXX = 1e12.
    rluxx = np.where(rlu < 9999.0, rlu / (hstar / 1e5 + f0), 1e12)

    # Soil: parallel SO2-like and O3-like pathways
    rgsx = 1.0 / (hstar / 1e5 / rgss + f0 / rgso)

    # Lower canopy: parallel pathways.
    # np.errstate: when rcls=rclo=inf the denominator is 0 → rclx=inf,
    # so 1/rclx=0, correctly excluding the pathway from the parallel sum.
    with np.errstate(divide="ignore", invalid="ignore"):
        rclx = 1.0 / (hstar / 1e5 / rcls + f0 / rclo)

    # ── Bulk surface resistance (Wesely Fig. 1 parallel network) ─────────
    r_surface = 1.0 / (1.0 / rixx + 1.0 / rluxx + 1.0 / (rac + rgsx) + 1.0 / (rdc + rclx))
    r_surface = np.clip(r_surface, 1.0, 9999.0)

    # ── Aerodynamic and quasi-laminar resistances (patch-independent) ─────
    ckustr = _VON_KARMAN * ustar
    reynolds = ustar * roughness_length / kin_visc
    corr = cz1 / obk
    z0_obk = roughness_length / obk

    if reynolds >= 0.1:
        # Aerodynamically rough surface
        if corr < 0.0:
            # Unstable
            d1 = (1.0 - 15.0 * corr) ** 0.5
            d2 = (1.0 - 15.0 * z0_obk) ** 0.5
            d3 = abs((d1 - 1.0) / (d1 + 1.0))
            d4 = abs((d2 - 1.0) / (d2 + 1.0))
            r_aero = (1.0 / ckustr) * np.log(d3 / d4)
        elif corr <= 1.0:
            # Weakly stable
            r_aero = (1.0 / ckustr) * (np.log(corr / z0_obk) + 5.0 * (corr - z0_obk))
        else:
            # Strongly stable
            r_aero = (1.0 / ckustr) * (5.0 * np.log(corr / z0_obk) + (corr - z0_obk))

        r_aero = float(np.clip(r_aero, 0.0, 1e4))

        r_quasi_laminar = (2.0 / ckustr) * (_D_AIR / diff_gas) ** 0.667
        c1x = r_aero + r_quasi_laminar + r_surface
    else:
        # Aerodynamically smooth surface
        r_aero = 1e4
        c1x = r_aero + r_surface  # no quasi-laminar term in original smooth branch

    # ── Sum contributions weighted by fractional cover ────────────────────
    vd = float(np.sum(0.001 * fractions * 1000.0 / c1x))
    return vd * 100.0  # convert m s⁻¹ → cm s⁻¹
