"""
Dust Emission Scheme: GOCART / Ginoux et al. (2001), WRF dust_opt=1
===================================================================
Target Parameterization for QUACS / CheMPAS-A

This module ports the WRF-Chem GOCART dust emission and driver used by
namelist option dust_opt=1. The flux loop, threshold
velocity, moisture cutoff, and bin parameters are pinned to the WRF source
below at commit 06d4240ae989cc3e50af412bb472df3d9048783c.

Fortran references (all line numbers at that commit):
    chem/module_gocart_dust.F
        subroutine source_du            -- per-bin flux loop + threshold
        subroutine gocart_dust_driver   -- gwet = smois / porosity(isltyp)
    phys/module_data_gocart_dust.F
        L2   ndust=5, ndcls=3, ndsrc=1
        L9   porosity (19 USGS soil types)
        L18  den_dust  = (2500,2650,2650,2650,2650)   [kg/m3]
        L19  reff_dust = (0.73,1.4,2.4,4.5,8.0)e-6    [m]
        L20  ipoint    = (3,2,2,2,2)                  [1-based class index]
        L21  frac_s    = (0.1,0.25,0.25,0.25,0.25)    [-] (sums to 1.1)
        L15  ch_dust(ndust,12) declared (no value here)
    chem/chemics_init.F
        L1694,1723,1742,1842,2029  ch_dust(:,:) = 0.8D-9
        (runtime value for the GOCART / dust_opt=1 aerosol-init CASEs)

WRF source_du per-bin flux, as a rate [kg m^-2 s^-1]:
    u_ts0_n = 0.13 * sqrt(rho_p_n * g * d_n / rho_a)
                   * sqrt(1 + 6e-7 / (rho_p_n * g * d_n^2.5))
                   / sqrt(1.928*(1331*(100*d_n)^1.56 + 0.38)^0.092 - 1)
    u_ts_n  = max(0, u_ts0_n * (1.2 + 0.2*log10(max(1e-3, gwet))))   if gwet < 0.5
              100 m/s otherwise
    F_p_n   = max(ch_dust[n,start_month] * frac_s_n * erod[ipoint_n]
                  * w10m^2 * (w10m - u_ts_n), 0)

Note there is NO lake-fraction term in dust_opt=1. The (1-frlake) factor is
GOCART2G (DustEmissionGOCART2G_revised), a different scheme; it is omitted here.

Migration notes (changes from the prior GOCART2G-flavoured version):
    1. u_t is now per-bin, shape [x, y, n_bins] (was a single [x, y] threshold
       collapsed to a bin mean). dust_emission() consumes threshold_velocity()
       output directly.
    2. frlake / lake_factor removed entirely (not in dust_opt=1).
    3. Erodibility is now the 3-class field erod[x, y, ndcls] with ipoint[n]
       selecting Sand/Silt/Clay per bin (was a single scalar S).
    4. Constants use WRF symbols: REFF_DUST, DEN_DUST, FRAC_S,
       IPOINT, CH_DUST (the latter is a [5,12] table).
    Code that imports this module must use the new dust_emission() signature
    and must not pass the removed frlake.

QUACS Contract
--------------
Standalone NumPy module. Runs as a script (python dust_ginoux.py) for the
box-model demo, or as a pytest suite (pytest dust_ginoux.py -v). The fast
vectorized dust_emission() is checked against a slow scalar transcription of
source_du to machine precision.

Inputs (dynamic, per time step from MPAS-A):
    w10m : np.ndarray [x, y]           10-m wind speed magnitude [m/s]
    gwet : np.ndarray [x, y]           volumetric soil moisture over porosity
                                       (saturation fraction) [-], NOT gravimetric
    u_t  : np.ndarray [x, y, n_bins]   per-bin dry threshold from threshold_velocity()

Inputs (static, loaded once at initialization):
    erod : np.ndarray [x, y, ndcls=3]  erodible fraction per class
                                       (index 0=Sand, 1=Silt, 2=Clay) [-], in [0,1]
    oro  : np.ndarray [x, y]           land mask (1.0 = land)
    rhoa : np.ndarray [x, y]           near-surface air density [kg/m^3]

Output:
    F_p : np.ndarray [x, y, n_bins]    surface emission flux [kg m^-2 s^-1]
                                       QUACS runtime: dx/dt = F_p / (rho_air * dz)
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field, fields

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.81  # m/s^2   (WRF driver g)
LAND = 1.0  # oro / xland value identifying land (xland < 1.5 in WRF)

# ---------------------------------------------------------------------------
# Bin and scheme parameters -- source: phys/module_data_gocart_dust.F
# ---------------------------------------------------------------------------
NDCLS = 3  # erodibility classes: 0=Sand, 1=Silt, 2=Clay  (L2: ndcls=3)

# Effective particle radius per bin [m]            (L19: reff_dust)
REFF_DUST = np.array([0.73, 1.4, 2.4, 4.5, 8.0]) * 1.0e-6  # [m]

# Particle bulk density per bin [kg/m^3]           (L18: den_dust)
# Bin 1 is 2500, bins 2-5 are 2650. This matches WRF dust_opt=1.
DEN_DUST = np.array([2500.0, 2650.0, 2650.0, 2650.0, 2650.0])  # [kg/m^3]

# Per-bin source mass fraction [-]                 (L21: frac_s)
# Sums to 1.1, not 1.0 -- this is the GOCART value, not a normalization error.
FRAC_S = np.array([0.1, 0.25, 0.25, 0.25, 0.25])  # sum = 1.1

# Bin -> erodibility class map, 0-based.           (L20: ipoint = 3,2,2,2,2)
# Fortran is 1-based (3,2,2,2,2): bin 1 -> Clay, bins 2-5 -> Silt; Sand unused.
IPOINT = np.array([3, 2, 2, 2, 2]) - 1  # = [2,1,1,1,1]

# Dust emission tuning coefficient [kg s^2 m^-5].
# The data module only DECLARES ch_dust(ndust,12) (L15); the value used on the
# dust_opt=1 path is assigned at runtime in chem/chemics_init.F:
#   ch_dust(:,:) = 0.8D-9  at L1694,1723,1742,1842,2029
#   (CASE GOCARTRACM_KPP / DUST / GOCART_SIMPLE / MOZCART_KPP aerosol init).
# NOTE: phys/module_dust_emis.F::bulk_dust_emis hardcodes 1.0e-9, but that is
# the QNIFA physics path, NOT gocart_dust_driver/source_du. Do not use it here.
# With CH_DUST in kg s2 m-5, F_p is returned directly in kg m-2 s-1.
CH_DUST = np.full((5, 12), 0.8e-9)  # ch_dust(bin, start_month)
CH_DUST_BULK_PHYS = 1.0e-9  # phys bulk_dust_emis only; not dust_opt=1

# Land mask cutoff. WRF sets the land flag with xland < 1.5
# (gocart_dust_driver: if(xland(i,j).lt.1.5) ... ilwi=1).
LAND_MAX = 1.5

# Soil porosity by USGS soil type [-]; source: phys/module_data_gocart_dust.F L9-12.
# gwet on the dust_opt=1 path is smois_top / porosity(isltyp), a saturation fraction.
POROSITY = np.array(
    [
        0.339,
        0.421,
        0.434,
        0.476,
        0.476,
        0.439,
        0.404,
        0.464,
        0.465,
        0.406,
        0.468,
        0.468,
        0.439,
        1.000,
        0.200,
        0.421,
        0.468,
        0.200,
        0.339,
    ]
)

LAYER_CUTOFF = 12.0
WET_CUTOFF, WET_THRESHOLD = 0.5, 100.0
MOISTURE_BASE, MOISTURE_LOG, MOISTURE_FLOOR = 1.2, 0.2, 1.0e-3
THRESHOLD_SCALE, THRESHOLD_COHESION, THRESHOLD_DRAG = 0.13, 6.0e-7, 1.928
THRESHOLD_SIZE, THRESHOLD_OFFSET = 1331.0, 0.38
THRESHOLD_SIZE_POWER, THRESHOLD_DRAG_POWER, THRESHOLD_DIAMETER_POWER = 1.56, 0.092, 2.5


@dataclass(frozen=True)
class GinouxSettings:
    reff_dust: tuple = field(
        default_factory=lambda: tuple(REFF_DUST),
        metadata={"source": "data.F:19", "what": "dust bin radii"},
    )
    den_dust: tuple = field(
        default_factory=lambda: tuple(DEN_DUST),
        metadata={"source": "data.F:18", "what": "dust density"},
    )
    frac_s: tuple = field(
        default_factory=lambda: tuple(FRAC_S),
        metadata={"source": "data.F:21", "what": "bin source fractions"},
    )
    ipoint: tuple = field(
        default_factory=lambda: tuple(IPOINT),
        metadata={"source": "data.F:20", "what": "bin soil class map"},
    )
    porosity: tuple = field(
        default_factory=lambda: tuple(POROSITY),
        metadata={"source": "data.F:9", "what": "soil porosity"},
    )
    ch_dust: tuple = field(
        default_factory=lambda: tuple(map(tuple, CH_DUST)),
        metadata={"source": "chemics_init.F:1694", "what": "monthly source coefficient"},
    )
    layer_cutoff: float = field(
        default=LAYER_CUTOFF,
        metadata={"source": "module_gocart_dust.F:87", "what": "thin layer wind cutoff"},
    )
    wet_cutoff: float = field(
        default=WET_CUTOFF,
        metadata={"source": "module_gocart_dust.F:192", "what": "wetness emission cutoff"},
    )
    wet_threshold: float = field(
        default=WET_THRESHOLD,
        metadata={"source": "module_gocart_dust.F:196", "what": "wet threshold wind"},
    )
    moisture_base: float = field(
        default=MOISTURE_BASE,
        metadata={"source": "module_gocart_dust.F:193", "what": "moisture threshold base"},
    )
    moisture_log: float = field(
        default=MOISTURE_LOG,
        metadata={"source": "module_gocart_dust.F:193", "what": "moisture logarithm scale"},
    )
    moisture_floor: float = field(
        default=MOISTURE_FLOOR,
        metadata={"source": "module_gocart_dust.F:193", "what": "wetness logarithm floor"},
    )
    threshold_scale: float = field(
        default=THRESHOLD_SCALE,
        metadata={"source": "module_gocart_dust.F:186", "what": "dry threshold scale"},
    )
    threshold_cohesion: float = field(
        default=THRESHOLD_COHESION,
        metadata={"source": "module_gocart_dust.F:186", "what": "dry threshold cohesion"},
    )
    threshold_drag: float = field(
        default=THRESHOLD_DRAG,
        metadata={"source": "module_gocart_dust.F:188", "what": "dry threshold drag"},
    )
    threshold_size: float = field(
        default=THRESHOLD_SIZE,
        metadata={"source": "module_gocart_dust.F:188", "what": "dry threshold size scale"},
    )
    threshold_offset: float = field(
        default=THRESHOLD_OFFSET,
        metadata={"source": "module_gocart_dust.F:188", "what": "dry threshold offset"},
    )
    threshold_size_power: float = field(
        default=THRESHOLD_SIZE_POWER,
        metadata={"source": "module_gocart_dust.F:188", "what": "dry threshold size exponent"},
    )
    threshold_drag_power: float = field(
        default=THRESHOLD_DRAG_POWER,
        metadata={"source": "module_gocart_dust.F:188", "what": "dry threshold drag exponent"},
    )
    threshold_diameter_power: float = field(
        default=THRESHOLD_DIAMETER_POWER,
        metadata={"source": "module_gocart_dust.F:187", "what": "dry threshold diameter exponent"},
    )

    def __post_init__(self):
        for name in ("reff_dust", "den_dust", "frac_s", "ipoint", "porosity", "ch_dust"):
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                value = value.copy()
                value.setflags(write=False)
                object.__setattr__(self, name, value)


def tunable_names(settings_class):
    return tuple(item.name for item in fields(settings_class))


def compute_gwet(smois_top, isltyp, settings=None):
    """
    Driver-level moisture variable used by source_du: volumetric soil moisture
    over porosity (a saturation fraction), NOT gravimetric moisture.

    Parameters
    ----------
    smois_top : np.ndarray   top-layer volumetric soil moisture [m^3/m^3]
    isltyp    : np.ndarray   integer USGS soil type, 1-based (1..19)

    Returns
    -------
    gwet : np.ndarray   smois_top / porosity(isltyp) [-]
    """

    settings = GinouxSettings() if settings is None else settings
    isltyp = np.asarray(isltyp)
    if np.any((isltyp < 1) | (isltyp > len(settings.porosity)) | (isltyp != np.floor(isltyp))):
        raise ValueError("[dust_ginoux] isltyp must be 1-based USGS soil type in [1, 19]")
    return smois_top / np.asarray(settings.porosity)[isltyp - 1]


# ---------------------------------------------------------------------------
# Threshold velocity: Marticorena & Bergametti (1995) / Bagnold (1941)
# Source: chem/module_gocart_dust.F, subroutine source_du (u_ts0).
# Verified SI transcription of the WRF CGS expression to ~1e-16.
# ---------------------------------------------------------------------------
def threshold_velocity(reff_dust, den_dust, rhoa, settings=None):
    """
    Dry-soil threshold wind speed per size bin and grid cell.

    Parameters
    ----------
    reff_dust : np.ndarray [n_bins]   effective particle radius per bin [m]
    den_dust  : np.ndarray [n_bins]   particle density per bin [kg/m^3]
    rhoa      : np.ndarray [x, y]     near-surface air density [kg/m^3]

    Returns
    -------
    u_t : np.ndarray [x, y, n_bins]   dry-soil threshold wind speed [m/s]
    """
    settings = GinouxSettings() if settings is None else settings
    d = (2.0 * np.asarray(reff_dust))[np.newaxis, np.newaxis, :]  # diameter [m], [1,1,n_bins]
    rp = np.asarray(den_dust)[np.newaxis, np.newaxis, :]  # [1,1,n_bins]
    ra = rhoa[:, :, np.newaxis]  # [x,y,1]

    u_t = (
        settings.threshold_scale
        * np.sqrt(rp * GRAV * d / ra)
        * np.sqrt(
            1.0 + settings.threshold_cohesion / (rp * GRAV * d**settings.threshold_diameter_power)
        )
        / np.sqrt(
            settings.threshold_drag
            * (
                settings.threshold_size * (100.0 * d) ** settings.threshold_size_power
                + settings.threshold_offset
            )
            ** settings.threshold_drag_power
            - 1.0
        )
    )
    return u_t  # [x, y, n_bins]


# ---------------------------------------------------------------------------
# Initialization: load static fields once at model start
# ---------------------------------------------------------------------------
def initialize(erod_file, oro_file, ndcls=NDCLS):
    """
    Load and validate static input fields. Called once before the time loop.

    Parameters
    ----------
    erod_file : str   Path to erod [x, y, ndcls], values in [0, 1].
    oro_file  : str   Path to land mask [x, y] (1.0 = land).

    Returns
    -------
    erod : np.ndarray [x, y, ndcls]
    oro  : np.ndarray [x, y]
    """
    for path, label in [(erod_file, "erod_file"), (oro_file, "oro_file")]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_ginoux] {label} not found: {path}")

    erod = np.load(erod_file)
    oro = np.load(oro_file)

    if erod.ndim != 3 or erod.shape[2] != ndcls:
        raise ValueError(f"[dust_ginoux] erod must be [x,y,{ndcls}], got {erod.shape}")
    if oro.shape != erod.shape[:2]:
        raise ValueError(f"[dust_ginoux] oro shape {oro.shape} != erod xy {erod.shape[:2]}")

    land_erod = erod[oro < LAND_MAX]
    if np.any((land_erod < 0) | (land_erod > 1) | ~np.isfinite(land_erod)):
        warnings.warn("[dust_ginoux] land erod outside [0,1]", UserWarning)

    return erod, oro


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# Source: chem/module_gocart_dust.F, subroutine source_du.
# ---------------------------------------------------------------------------
def dust_emission(w10m, u_t, gwet, oro, erod, month=1, settings=None):
    """
    Return source_du flux rate [kg m-2 s-1] on land; bems equals flux * dt.

    ch_dust has shape [5,12], indexed by bin and month; [x,y,5,12] allows
    paper samples to vary C by cell. WRF selects the run's start_month once
    and uses that column at every time step. frac_s and ipoint are the fixed
    WRF constants FRAC_S and IPOINT (data module L20-21), not inputs.
    """
    settings = GinouxSettings() if settings is None else settings
    ch_dust = settings.ch_dust
    shape = np.shape(oro)
    if (
        np.shape(w10m) != shape
        or np.shape(gwet) != shape
        or np.shape(erod) != (*shape, 3)
        or np.shape(u_t) != (*shape, 5)
    ):
        raise ValueError("[dust_ginoux] inconsistent grid or bin shapes")
    if np.shape(ch_dust) not in ((5, 12), (*shape, 5, 12)) or not 1 <= month <= 12:
        raise ValueError("[dust_ginoux] ch_dust must be [5,12] or [x,y,5,12], month in 1..12")
    flux = np.zeros((*shape, 5))
    land = np.asarray(oro) < LAND_MAX
    if not land.any():
        return flux
    wet = np.asarray(gwet)[land]
    erosion = np.asarray(erod)[land]
    if np.any((erosion < 0) | (erosion > 1) | ~np.isfinite(erosion)):
        warnings.warn("[dust_ginoux] land erod outside [0,1]", UserWarning)
    if np.any((wet < 0) | (wet > 1) | ~np.isfinite(wet)):
        warnings.warn("[dust_ginoux] land gwet outside [0,1]", UserWarning)
    # fmax, not maximum: gfortran MAX(0, NaN) returns 0 (F L193), fmax gives the non-NaN argument.
    threshold = np.where(
        wet[:, None] < settings.wet_cutoff,
        np.fmax(
            np.asarray(u_t)[land]
            * (
                settings.moisture_base
                + settings.moisture_log * np.log10(np.fmax(wet[:, None], settings.moisture_floor))
            ),
            0.0,
        ),
        settings.wet_threshold,
    )
    wind = np.asarray(w10m)[land, None]
    table = np.asarray(ch_dust)
    coefficient = table[:, month - 1] if table.ndim == 2 else table[land][:, :, month - 1]
    if np.any((coefficient < 0) | ~np.isfinite(coefficient)):
        warnings.warn("[dust_ginoux] land ch_dust is negative or nonfinite", UserWarning)
    if np.any((wind < 0) | ~np.isfinite(wind)):
        warnings.warn("[dust_ginoux] land wind is negative or nonfinite", UserWarning)
    raw = (
        coefficient
        * np.asarray(settings.frac_s)
        * erosion[:, settings.ipoint]
        * wind**2
        * (wind - threshold)
    )
    flux[land] = np.maximum(raw, 0.0)
    return flux


def ginoux_dust_driver(
    dust,
    rho_phy,
    dz8w,
    dt,
    u10,
    v10,
    u_phy,
    v_phy,
    smois,
    isltyp,
    oro,
    erod,
    emis_dust,
    dustin,
    chem_opt=0,
    month=1,
    settings=None,
):
    """Port gocart_dust_driver for five dust tracers in ug/kg.

    Fields use WRF layout: dust [x,level,y,5], rho_phy/dz8w/u_phy/v_phy
    [x,level,y], smois [x,soil,y], surface fields [x,y]. Return flux in
    kg m-2 s-1 and updated dust, emis_dust and dustin. WRF selects the run's
    start month. Sea cells retain their starting values.
    """
    settings = GinouxSettings() if settings is None else settings
    land = np.asarray(oro) < LAND_MAX
    updated = np.array(dust, dtype=float, copy=True)
    emitted = np.array(emis_dust, dtype=float, copy=True)
    incoming = np.array(dustin, dtype=float, copy=True)
    flux = np.zeros((*np.shape(oro), 5))
    if not land.any():
        return dict(flux=flux, dust=updated, emis_dust=emitted, dustin=incoming)
    soil = np.asarray(isltyp)[land]
    if np.any((soil < 1) | (soil > 19) | (soil != np.floor(soil))):
        raise ValueError("[dust_ginoux] land isltyp must be integer 1..19")
    wind = np.sqrt(np.asarray(u10)[land] ** 2 + np.asarray(v10)[land] ** 2)
    thin = np.asarray(dz8w)[:, 0, :][land] < settings.layer_cutoff
    model_wind = np.sqrt(
        np.asarray(u_phy)[:, 0, :][land] ** 2 + np.asarray(v_phy)[:, 0, :][land] ** 2
    )
    wind = np.where(thin, model_wind, wind)
    wet = np.asarray(smois)[:, 0, :][land] / np.asarray(settings.porosity)[soil.astype(int) - 1]
    density = np.asarray(rho_phy)[:, 0, :][land]
    threshold = threshold_velocity(
        settings.reff_dust, settings.den_dust, density[:, None], settings
    )
    flux[land] = dust_emission(
        wind[:, None],
        threshold,
        wet[:, None],
        np.ones((len(wind), 1)),
        np.asarray(erod)[land, None, :],
        month=month,
        settings=settings,
    )[:, 0, :]
    mass = flux[land] * dt
    current = updated[:, 0, :, :]
    start = np.full_like(current[land], 1e-16) if chem_opt in (2, 11) else current[land]
    final = start + mass / np.asarray(dz8w)[:, 0, :][land, None] / density[:, None] * 1e9
    if chem_opt in (2, 11):
        incoming[land] = final
    else:
        current[land] = final
    emitted[land] += mass
    return dict(flux=flux, dust=updated, emis_dust=emitted, dustin=incoming)


# ===========================================================================
# Box-model demo  --  run with: python dust_ginoux.py
# ===========================================================================
