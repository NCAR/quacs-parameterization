"""
Dust Emission Scheme: AFWA / GOCART-AFWA (WRF-Chem dust_opt=3)
==============================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    LeGrand, S. L., Polashenski, C., Letcher, T. W., Creighton, G. A.,
    Peckham, S. E., and Cetola, J. D. (2019). The AFWA dust emission scheme
    for the GOCART aerosol model in WRF-Chem v3.8.1.
    Geosci. Model Dev., 12, 131-166. https://doi.org/10.5194/gmd-12-131-2019

    Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust
    cycle: 1. Design of a soil-derived dust emission scheme.
    J. Geophys. Res., 100, 16415-16430. https://doi.org/10.1029/95JD00690
    (hereon MB95: saltation flux Eq. 28, surface-area weighting Eq. 32,
     sandblasting efficiency beta.)

    Iversen, J. D. and White, B. R. (1982). Saltation threshold on Earth,
    Mars and Venus. Sedimentology, 29, 111-119. (dry threshold form)

    Fecan, F., Marticorena, B., Bergametti, G. (1999). Parametrization of the
    increase of the aeolian erosion threshold wind friction velocity due to
    soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157.
    (moisture correction to threshold friction velocity.)

    Kok, J. F. (2011). A scaling theory for the size distribution of emitted
    dust aerosols. PNAS, 108, 1016-1021. (emitted-dust bin distribution.)

    Fortran reference: WRF-Chem (wrf-model/WRF, master)
      chem/module_gocart_dust_afwa.F   subroutine source_dust
      phys/module_data_gocart_dust.F   all bin/soil parameters
      Registry/registry.chem           dust_alpha, dust_gamma, smtune, ustune
    Scheme author credit (Fortran): Sandra Jones (AER/AFWA) and
    Glenn Creighton (AFWA).

QUACS Contract
--------------
Standalone Python module. No MPAS-A, no config files, no external
dependencies beyond NumPy. Runs as a script (python dust_afwa.py)
or as a pytest test suite (pytest dust_afwa.py -v).

Pre-API note:
    Once the QUACS scripting API is finalized, the expected signature is:
        @scheme(category='dust_emission')
        def dust_emission_afwa(state: ModelState, config: Config) -> SurfaceFlux:
    Until then, plain NumPy arrays are used. This script is the scientific
    specification that drives that API design.

How AFWA differs from Ginoux/GOCART (dust_opt=1, see dust_ginoux.py)
-------------------------------------------------------------------
Ginoux is an empirical wind-power law driven by 10-m wind speed (w10m) and a
single moisture-adjusted threshold velocity. AFWA is a physically based
saltation/sandblasting scheme:
    1. It is driven by FRICTION VELOCITY u* (not w10m).
    2. Horizontal saltation flux is built up over 9 saltation (sand) bins
       using the MB95 white-noise saltation law, weighted by a soil-class
       surface-area distribution.
    3. The horizontal-to-vertical conversion uses a clay-dependent
       sandblasting efficiency beta and a global tuning constant alpha.
    4. The emitted mass is split into 5 dust bins by a fixed Kok (2011)
       lognormal distribution (distr_dust), not by a per-bin source fraction.

Physics summary (source_dust):
    Per saltation bin n (smx = 9):
        dry threshold (Iversen-White, CGS):
            u_ts0_n = 0.13e-2 * sqrt(den_n*g*diam_n/rhoa)
                      * sqrt(1 + 0.006/(den_n*g*diam_n^2.5))
                      / sqrt(1.928*(1331*diam_n^1.56 + 0.38)^0.092 - 1)
        moisture correction (gravimetric, Fecan-type):
            u_ts_n = u_ts0_n * sqrt(1 + 1.21*(gravsm - drylimit)^0.68)   if gravsm > drylimit
                   = u_ts0_n                                              otherwise
        saltation flux (MB95), only where u* > u_ts_n:
            salt_n = cmb * ds_rel_n * (airden/g0) * u*^3
                     * (1 + u_ts_n/u*) * (1 - u_ts_n^2/u*^2)
    Total vertical (emitted) mass flux:
        beta  = min(10^(13.6*clay - 6.0), betamax)
        emit  = (sum_n salt_n) * erod^gamma * alpha * beta            [kg m^-2 s^-1]
    Split into 5 dust bins:
        F_p[m] = emit * distr_dust[m]                                 [kg m^-2 s^-1]

    Emission is zeroed over water (oro >= 1.5), where erodibility is 0, where
    roughness znt > 0.2 m, where snow depth > 0.01 m, and over non-erodible
    soil types (bedrock/lava/land-ice).

Inputs (dynamic, per time step from MPAS-A):
    ustar  : np.ndarray [x, y]     surface friction velocity u* [m/s]
    volsm  : np.ndarray [x, y]     volumetric soil moisture [m^3/m^3]
    airden : np.ndarray [x, y]     near-surface air density [kg/m^3]
    znt    : np.ndarray [x, y]     surface roughness length [m]
    snowh  : np.ndarray [x, y]     snow depth (water-equivalent height) [m]

Inputs (static, loaded once at initialization):
    erod   : np.ndarray [x, y, 3]  Ginoux erodibility bins; summed internally
    clay   : np.ndarray [x, y]     soil clay mass fraction [-], in [0, 1]
    sand   : np.ndarray [x, y]     soil sand mass fraction [-], in [0, 1]
    oro    : np.ndarray [x, y]     land-water mask (WRF xland < 1.5 = land)
    vegfra : np.ndarray [x, y]     vegetation fraction [%], for dust_veg=1
    lai_vegmask: np.ndarray [x, y] MODIS vegetation mask, for dust_veg=2
    isltyp : np.ndarray [x, y]     soil type index (required; excludes 15/16/18)

Configurable (WRF Registry namelist defaults, registry.chem):
    dust_veg : int   vegetation option (default 0: Ginoux; 1: greenfrac; 2: LAI)
    alpha  : float   global tuning constant [m^-1]   (dust_alpha,  default 1.0)
    gamma  : float   erodibility exponent  [-]        (dust_gamma,  default 1.0)
    smtune : float   soil-moisture tuning  [-]        (dust_smtune, default 1.0)
    ustune : float   u* tuning             [-]        (dust_ustune, default 1.0)
    dust_smois, dust_dsr, dust_soils : int   WRF selectors, default 0
    sf_surface_physics : int   WRF LSM number (2 Noah, 3 RUC, 7 PX)
    cmb is fixed at 1.0 (afwa.F L488); no caller override.

Output:
    F_p : np.ndarray [x, y, n_dust]   surface emission flux [kg m^-2 s^-1]
                                      QUACS runtime converts to tendency:
                                      dx/dt = F_p / (airden * dz)
                                      (matches the AFWA Fortran update
                                       tc += dsrc / dz_lowest / airden).

Note on the erod**gamma + alpha + beta chain:
    alpha [m^-1] is the dimensional bridge that turns the horizontal saltation
    flux [kg m^-1 s^-1] into a vertical emission flux [kg m^-2 s^-1]. gamma
    reshapes the spatial weighting by erodibility (1.0 = neutral). beta is the
    physical sandblasting efficiency, rising with clay content and capped at
    betamax. All three are bin-independent and factor out of the saltation sum.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field, replace
from dataclasses import fields as dataclass_fields

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.81  # m/s^2; WRF v4.8.0 share/module_model_constants.F:17; CGS g = GRAV*100
LAND = 1.0  # oro value identifying land cells
ZNT_MAX = 0.2  # m   max roughness length that may loft (afwa.F kludge)
SNOW_MAX = 0.01  # m   snow depth above which lofting is suppressed
NONERODIBLE_SOIL = (15, 16, 18)  # bedrock, land-ice, lava (isltyp indices)

# ---------------------------------------------------------------------------
# Bin and soil parameters - source: phys/module_data_gocart_dust.F
# ---------------------------------------------------------------------------

# Number of saltation (sand) bins and emitted-dust bins
NGSALT = 9  # ngsalt
NDUST = 5  # ndust

# Saltation reference radii [m]  (reff_salt)
REFF_SALT = np.array([0.71, 1.37, 2.63, 5.00, 9.50, 18.1, 34.5, 65.5, 125.0]) * 1.0e-6  # [m]

# Saltation particle density per bin [kg/m^3]  (den_salt)
DEN_SALT = np.array(
    [2500.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0]
)  # [kg/m^3]

# Soil-class pointer per saltation bin (spoint): 1=clay, 2=silt, 3=sand
SPOINT = np.array([1, 2, 2, 2, 2, 2, 3, 3, 3])

# Soil-class mass fraction per saltation bin (frac_salt)
FRAC_SALT = np.array([1.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.333, 0.333, 0.333])

# Emitted-dust bin distribution (distr_dust) - hardcoded Kok (2011) lognormal.
# Source comment in afwa.F: precomputed and frozen rather than recomputed each
# step. Sums to 1.0000.
DISTR_DUST = np.array([1.074e-1, 1.012e-1, 2.078e-1, 4.817e-1, 1.019e-1])

# Emitted-dust reference radii [m]  (reff_dust) - carried for documentation
REFF_DUST = np.array([0.73, 1.4, 2.4, 4.5, 8.0]) * 1.0e-6  # [m]

# Soil porosity by WRF soil-type index (porosity, 19 classes).
# Used to convert volumetric -> gravimetric soil moisture.
POROSITY_TABLE = np.array(
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

# ---------------------------------------------------------------------------
# Tuning / efficiency constants - source: module_gocart_dust_afwa.F + Registry
# ---------------------------------------------------------------------------
BETAMAX = 5.25e-4  # max sandblasting mass efficiency [-]
CMB = 1.0  # saltation proportionality (Marticorena 1997).
# White (1979) value 2.61 is present-but-commented in afwa.F.

ALPHA_DEFAULT = 1.0  # dust_alpha  [m^-1]  (registry.chem)
GAMMA_DEFAULT = 1.0  # dust_gamma  [-]
SMTUNE_DEFAULT = 1.0  # dust_smtune [-]
USTUNE_DEFAULT = 1.0  # dust_ustune [-]  (recommend 0.7 for PX LSM)

VEG_CUTOFF = 5.0
BETA_CLAY, BETA_OFFSET = 13.6, -6.0
DRY_RUC_SCALE, DRY_RUC_CLAY, DRY_RUC_OFFSET, DRY_RUC_POWER = 0.035, 13.52, 3.53, 2.68
DRY_NOAH_SCALE, DRY_NOAH_CLAY, DRY_NOAH_OFFSET, DRY_NOAH_POWER = 0.0756, 15.127, 3.09, 2.3211
DRY_BASE_QUAD, DRY_BASE_LINEAR = 14.0, 17.0
SOIL_DENS_SAND, SOIL_DENS_CLAY = 2.65, 2.50
MOIST_SCALE, MOIST_POWER = 1.21, 0.68
IW_SCALE, IW_COHESION, IW_DRAG = 0.13e-2, 0.006, 1.928
IW_SIZE, IW_OFFSET, IW_SIZE_POWER, IW_DRAG_POWER, IW_DIAM_POWER = 1331.0, 0.38, 1.56, 0.092, 2.5
LOFT_HEIGHT = 10.0
SURFACE_GEOMETRY = 0.75
EXTINCTION = (1.470e-6, 7.877e-7, 4.623e-7, 2.429e-7, 1.387e-7)
VISIBILITY_SCALE = 3.912


@dataclass(frozen=True)
class AfwaSettings:
    alpha: float = field(
        default=ALPHA_DEFAULT,
        metadata={"source": "registry.chem:dust_alpha", "what": "emission scale"},
    )
    gamma: float = field(
        default=GAMMA_DEFAULT,
        metadata={"source": "registry.chem:dust_gamma", "what": "erodibility exponent"},
    )
    smtune: float = field(
        default=SMTUNE_DEFAULT,
        metadata={"source": "registry.chem:dust_smtune", "what": "soil moisture scale"},
    )
    ustune: float = field(
        default=USTUNE_DEFAULT,
        metadata={"source": "registry.chem:dust_ustune", "what": "friction velocity scale"},
    )
    znt_max: float = field(
        default=ZNT_MAX, metadata={"source": "afwa.F:232", "what": "roughness cutoff"}
    )
    snow_max: float = field(
        default=SNOW_MAX, metadata={"source": "afwa.F:246", "what": "snow cutoff"}
    )
    nonerodible_soil: tuple = field(
        default=NONERODIBLE_SOIL, metadata={"source": "afwa.F:238", "what": "masked soil types"}
    )
    reff_salt: tuple = field(
        default_factory=lambda: tuple(REFF_SALT),
        metadata={"source": "data.F:23", "what": "saltation radii"},
    )
    den_salt: tuple = field(
        default_factory=lambda: tuple(DEN_SALT),
        metadata={"source": "data.F:24", "what": "saltation densities"},
    )
    spoint: tuple = field(
        default_factory=lambda: tuple(SPOINT),
        metadata={"source": "data.F:25", "what": "soil class per bin"},
    )
    frac_salt: tuple = field(
        default_factory=lambda: tuple(FRAC_SALT),
        metadata={"source": "data.F:26", "what": "soil fraction per bin"},
    )
    distr_dust: tuple = field(
        default_factory=lambda: tuple(DISTR_DUST),
        metadata={"source": "data.F:28", "what": "emitted dust split"},
    )
    porosity_table: tuple = field(
        default_factory=lambda: tuple(POROSITY_TABLE),
        metadata={"source": "data.F:9", "what": "soil porosity"},
    )
    betamax: float = field(
        default=BETAMAX, metadata={"source": "afwa.F:473", "what": "efficiency cap"}
    )
    cmb: float = field(default=CMB, metadata={"source": "afwa.F:488", "what": "saltation scale"})
    beta_clay: float = field(
        default=BETA_CLAY, metadata={"source": "afwa.F:595", "what": "clay efficiency slope"}
    )
    beta_offset: float = field(
        default=BETA_OFFSET, metadata={"source": "afwa.F:595", "what": "clay efficiency offset"}
    )
    veg_cutoff: float = field(
        default=VEG_CUTOFF, metadata={"source": "afwa.F:172", "what": "vegetation mask cutoff"}
    )
    dry_ruc_scale: float = field(
        default=DRY_RUC_SCALE, metadata={"source": "afwa.F:268-278", "what": "RUC dry limit scale"}
    )
    dry_ruc_clay: float = field(
        default=DRY_RUC_CLAY,
        metadata={"source": "afwa.F:268-278", "what": "RUC dry limit clay slope"},
    )
    dry_ruc_offset: float = field(
        default=DRY_RUC_OFFSET,
        metadata={"source": "afwa.F:268-278", "what": "RUC dry limit offset"},
    )
    dry_ruc_power: float = field(
        default=DRY_RUC_POWER,
        metadata={"source": "afwa.F:268-278", "what": "RUC dry limit exponent"},
    )
    dry_noah_scale: float = field(
        default=DRY_NOAH_SCALE,
        metadata={"source": "afwa.F:268-278", "what": "Noah dry limit scale"},
    )
    dry_noah_clay: float = field(
        default=DRY_NOAH_CLAY,
        metadata={"source": "afwa.F:268-278", "what": "Noah dry limit clay slope"},
    )
    dry_noah_offset: float = field(
        default=DRY_NOAH_OFFSET,
        metadata={"source": "afwa.F:268-278", "what": "Noah dry limit offset"},
    )
    dry_noah_power: float = field(
        default=DRY_NOAH_POWER,
        metadata={"source": "afwa.F:268-278", "what": "Noah dry limit exponent"},
    )
    dry_base_quad: float = field(
        default=DRY_BASE_QUAD,
        metadata={"source": "afwa.F:268-278", "what": "base dry limit quadratic"},
    )
    dry_base_linear: float = field(
        default=DRY_BASE_LINEAR,
        metadata={"source": "afwa.F:268-278", "what": "base dry limit linear"},
    )
    soil_dens_sand: float = field(
        default=SOIL_DENS_SAND, metadata={"source": "afwa.F:257", "what": "sand density conversion"}
    )
    soil_dens_clay: float = field(
        default=SOIL_DENS_CLAY, metadata={"source": "afwa.F:257", "what": "clay density conversion"}
    )
    moist_scale: float = field(
        default=MOIST_SCALE, metadata={"source": "afwa.F:560", "what": "wet threshold scale"}
    )
    moist_power: float = field(
        default=MOIST_POWER, metadata={"source": "afwa.F:566", "what": "wet threshold exponent"}
    )
    iw_scale: float = field(
        default=IW_SCALE, metadata={"source": "afwa.F:542-544", "what": "dry threshold scale"}
    )
    iw_cohesion: float = field(
        default=IW_COHESION, metadata={"source": "afwa.F:542-544", "what": "dry threshold cohesion"}
    )
    iw_drag: float = field(
        default=IW_DRAG, metadata={"source": "afwa.F:542-544", "what": "dry threshold drag"}
    )
    iw_size: float = field(
        default=IW_SIZE, metadata={"source": "afwa.F:542-544", "what": "dry threshold size scale"}
    )
    iw_offset: float = field(
        default=IW_OFFSET, metadata={"source": "afwa.F:542-544", "what": "dry threshold offset"}
    )
    iw_size_power: float = field(
        default=IW_SIZE_POWER,
        metadata={"source": "afwa.F:542-544", "what": "dry threshold size exponent"},
    )
    iw_drag_power: float = field(
        default=IW_DRAG_POWER,
        metadata={"source": "afwa.F:542-544", "what": "dry threshold drag exponent"},
    )
    iw_diam_power: float = field(
        default=IW_DIAM_POWER,
        metadata={"source": "afwa.F:542-544", "what": "dry threshold diameter exponent"},
    )
    loft_height: float = field(
        default=LOFT_HEIGHT, metadata={"source": "afwa.F:327", "what": "lofting reference height"}
    )
    extinction: tuple = field(
        default=EXTINCTION, metadata={"source": "afwa.F:356-357", "what": "bin extinction factors"}
    )
    visibility_scale: float = field(
        default=VISIBILITY_SCALE,
        metadata={"source": "afwa.F:356-357", "what": "visibility numerator"},
    )
    surface_geometry: float = field(
        default=SURFACE_GEOMETRY,
        metadata={
            "source": "afwa.F:508",
            "what": "sphere geometry factor 3/4 for area per unit mass; default from sphere geometry",
        },
    )

    def __post_init__(self):
        for name in (
            "reff_salt",
            "den_salt",
            "spoint",
            "frac_salt",
            "distr_dust",
            "porosity_table",
            "extinction",
        ):
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                value = value.copy()
                value.setflags(write=False)
                object.__setattr__(self, name, value)


def tunable_names(settings_class):
    return tuple(item.name for item in dataclass_fields(settings_class))


# ---------------------------------------------------------------------------
# Dry threshold friction velocity: Iversen & White (1982) form, CGS internally
# ---------------------------------------------------------------------------
def threshold_friction_velocity(reff_salt, den_salt, airden, settings=None):
    """
    Dry-soil threshold friction velocity u_ts0 per saltation bin and cell.

    Faithful to module_gocart_dust_afwa.F (source_dust), which works in CGS:
        den  = den_salt * 1e-3   (g cm^-3)
        diam = 2*reff_salt * 1e2 (cm)
        rhoa = airden * 1e-3     (g cm^-3)
        g    = g0 * 1e2          (cm s^-2)
        u_ts0 = 0.13e-2 * sqrt(den*g*diam/rhoa)
                * sqrt(1 + 0.006/(den*g*diam^2.5))
                / sqrt(1.928*(1331*diam^1.56 + 0.38)^0.092 - 1)
    The leading 1e-2 returns u_ts0 in m/s.

    Unlike the fine-only GOCART bins, AFWA's 9 saltation bins span
    0.71-125 um and cross the Iversen-White minimum (~30-100 um), so u_ts0
    is U-shaped across bins: it decreases, reaches a minimum, then rises.

    Parameters
    ----------
    reff_salt : np.ndarray [smx]    saltation reference radius [m]
    den_salt  : np.ndarray [smx]    saltation particle density [kg/m^3]
    airden    : np.ndarray [x, y]   near-surface air density [kg/m^3]

    Returns
    -------
    u_ts0 : np.ndarray [x, y, smx]   dry threshold friction velocity [m/s]
    """
    settings = AfwaSettings() if settings is None else settings
    den = (np.asarray(den_salt) * 1.0e-3)[np.newaxis, np.newaxis, :]  # [1,1,smx] g/cm^3
    diam = (2.0 * np.asarray(reff_salt) * 1.0e2)[np.newaxis, np.newaxis, :]  # [1,1,smx] cm
    rhoa = (airden * 1.0e-3)[:, :, np.newaxis]  # [x,y,1] g/cm^3
    g = GRAV * 1.0e2  # cm/s^2

    u_ts0 = (
        settings.iw_scale
        * np.sqrt(den * g * diam / rhoa)
        * np.sqrt(1.0 + settings.iw_cohesion / (den * g * diam**settings.iw_diam_power))
        / np.sqrt(
            settings.iw_drag
            * (settings.iw_size * diam**settings.iw_size_power + settings.iw_offset)
            ** settings.iw_drag_power
            - 1.0
        )
    )
    return u_ts0  # [x, y, smx]


# ---------------------------------------------------------------------------
# Saltation surface-area weighting: MB95 Eq. 32
# ---------------------------------------------------------------------------
def saltation_area_fraction(clay, sand, settings=None):
    """
    Relative saltation surface-area distribution ds_rel per bin (MB95 Eq. 32).

    Faithful to source_dust:
        dmass_n    = massfrac(spoint_n) * frac_salt_n
        dsurface_n = 0.75 * dmass_n / (den_salt_n * reff_salt_n)
        ds_rel_n   = dsurface_n / sum_n dsurface_n
    where massfrac = (clay, silt, sand) and silt = 1 - (clay + sand).
    Nonnegative fractions with clay + sand <= 1 always have positive area,
    including pure silt (clay = sand = 0). Like Fortran L519-522, division
    is unguarded: invalid fractions can yield nonfinite or unphysical values.

    Parameters
    ----------
    clay : np.ndarray [x, y]   clay mass fraction [-]
    sand : np.ndarray [x, y]   sand mass fraction [-]

    Returns
    -------
    ds_rel : np.ndarray [x, y, smx]   relative surface-area fraction [-]
    """
    settings = AfwaSettings() if settings is None else settings
    silt = 1.0 - (clay + sand)
    massfrac = np.stack([clay, silt, sand], axis=-1)  # [x,y,3] 0=clay,1=silt,2=sand
    mf_bin = massfrac[:, :, np.asarray(settings.spoint) - 1]  # [x,y,smx] gather by soil class

    dmass = mf_bin * np.asarray(settings.frac_salt)[np.newaxis, np.newaxis, :]
    dsurface = (
        settings.surface_geometry
        * dmass
        / (np.asarray(settings.den_salt) * np.asarray(settings.reff_salt))[
            np.newaxis, np.newaxis, :
        ]
    )
    stotal = dsurface.sum(axis=-1, keepdims=True)

    ds_rel = dsurface / stotal  # Fortran L519-522: no guard for invalid soil
    return ds_rel  # [x, y, smx]


# ---------------------------------------------------------------------------
# Initialization: load static external datasets once at model start
# ---------------------------------------------------------------------------
def initialize(
    erod_file,
    clay_file,
    sand_file,
    isltyp_file,
    oro_file,
    dust_veg=0,
    vegfra_file=None,
    lai_vegmask_file=None,
    dust_dsr=0,
    erod_dri_file=None,
    dust_soils=0,
    clay_nga_file=None,
    sand_nga_file=None,
    settings=None,
):
    """
    Load static fields and check file existence, array shapes and land soil types.

    Returns soil type and raw erodibility; it does not apply
    soil or vegetation exclusions. The timestep function constructs the masks.
    Soil types are checked only where xland < 1.5; water values are retained.
    Erodibility and texture values are loaded without new range restrictions.
    dust_veg other than 1 or 2 needs no vegetation file; option 1 requires vegfra_file (percent),
    and option 2 requires lai_vegmask_file. dust_dsr=1 requires erod_dri_file;
    dust_soils=1 requires clay_nga_file and sand_nga_file. Only selected files
    are loaded; negative regional sentinels are retained for per-cell fallback.

    Parameters
    ----------
    erod_file : str   Path to raw Ginoux erodibility [x, y, 3].
    clay_file, sand_file : str   Paths to soil fractions [x, y].
    isltyp_file : str   Path to 1-based WRF soil types [x, y].
    oro_file : str   Path to WRF xland [x, y], land where xland < 1.5.
    dust_veg : int   WRF vegetation selector; 1 greenfrac, 2 LAI, else Ginoux.
    vegfra_file, lai_vegmask_file : str   Files required only by their option.
    dust_dsr, dust_soils : int   Select DRI/NGA datasets only when equal to 1.
    erod_dri_file, clay_nga_file, sand_nga_file : str   Selected regional files.

    Returns
    -------
    erod, clay, sand, isltyp, oro, vegfra, lai_vegmask, erod_dri, clay_nga, sand_nga : tuple
        erod and erod_dri have shape [x, y, 3]; other arrays have shape [x, y].
        Unselected datasets are None. Pass the chosen options to the timestep.
    """
    settings = AfwaSettings() if settings is None else settings
    files = [erod_file, clay_file, sand_file, isltyp_file, oro_file]
    labels = ["erod_file", "clay_file", "sand_file", "isltyp_file", "oro_file"]
    for path, label in zip(files, labels):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_afwa] {label} not found: {path}")

    erod = np.load(erod_file)
    clay = np.load(clay_file)
    sand = np.load(sand_file)
    isltyp = np.load(isltyp_file)
    oro = np.load(oro_file)

    if erod.ndim != 3 or erod.shape[-1] != 3:
        raise ValueError("[dust_afwa] erod must have shape [x, y, 3]")
    ref = erod.shape[:2]
    for arr, name in [(clay, "clay"), (sand, "sand"), (isltyp, "isltyp"), (oro, "oro")]:
        if arr.shape != ref:
            raise ValueError(f"[dust_afwa] {name} shape {arr.shape} != erod grid shape {ref}")

    land = oro < 1.5
    soil = isltyp[land]
    if ((soil < 1) | (soil > len(settings.porosity_table)) | (soil != np.floor(soil))).any():
        raise ValueError("[dust_afwa] land isltyp must be an integer WRF soil type from 1 to 19")
    _warn_land_inputs(erod[land], clay[land], sand[land])
    needed = []
    if dust_veg == 1:
        needed.append(("vegfra", vegfra_file, ref))
    elif dust_veg == 2:
        needed.append(("lai_vegmask", lai_vegmask_file, ref))
    if dust_dsr == 1:
        needed.append(("erod_dri", erod_dri_file, erod.shape))
    if dust_soils == 1:
        needed.extend([("clay_nga", clay_nga_file, ref), ("sand_nga", sand_nga_file, ref)])
    loaded = {}
    for label, path, shape in needed:
        if path is None:
            raise ValueError(f"[dust_afwa] selected option requires {label}_file")
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_afwa] {label}_file not found: {path}")
        arr = np.load(path)
        if arr.shape != shape:
            raise ValueError(f"[dust_afwa] {label} shape {arr.shape} != expected {shape}")
        loaded[label] = arr

    return (
        erod,
        clay,
        sand,
        isltyp,
        oro,
        loaded.get("vegfra"),
        loaded.get("lai_vegmask"),
        loaded.get("erod_dri"),
        loaded.get("clay_nga"),
        loaded.get("sand_nga"),
    )


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def _afwa_land_fields(fields, land):
    """Pack only land cells into the existing two-dimensional flux interface."""
    packed = {
        name: value
        if value is None or np.ndim(value) == 0
        else np.asarray(value)[land][:, None, ...]
        for name, value in fields.items()
        if name != "settings"
    }
    if "settings" in fields:
        settings = fields["settings"]
        spatial = {
            item.name: np.asarray(getattr(settings, item.name))[land][:, None]
            for item in dataclass_fields(settings)
            if np.shape(getattr(settings, item.name)) == land.shape
        }
        packed["settings"] = replace(settings, **spatial) if spatial else settings
    return packed


def _warn_land_inputs(erod, clay, sand, volsm=None):
    """Warn on impossible land inputs while retaining WRF's arithmetic."""
    if np.any((erod < 0) | (erod > 1) | ~np.isfinite(erod)):
        warnings.warn("[dust_afwa] land erod outside [0,1]", UserWarning)
    if np.any(
        (clay < 0)
        | (clay > 1)
        | ~np.isfinite(clay)
        | (sand < 0)
        | (sand > 1)
        | ~np.isfinite(sand)
        | (clay + sand > 1)
    ):
        warnings.warn("[dust_afwa] land clay or sand outside soil fractions", UserWarning)
    if volsm is not None and np.any((volsm < 0) | (volsm > 1) | ~np.isfinite(volsm)):
        warnings.warn("[dust_afwa] land volsm outside [0,1]", UserWarning)


def _afwa_surface_fields(
    erod,
    clay,
    sand,
    dust_dsr=0,
    erod_dri=None,
    dust_soils=0,
    clay_nga=None,
    sand_nga=None,
    dust_veg=0,
    vegfra=None,
    lai_vegmask=None,
    settings=None,
):
    """WRF's regional selection and vegetation mask (F L154-223)."""
    settings = AfwaSettings() if settings is None else settings
    erod_total = erod.sum(axis=-1)
    if dust_dsr == 1:
        if erod_dri is None:
            raise ValueError("[dust_afwa] dust_dsr=1 requires erod_dri")
        erod_total = np.where(erod_dri[..., 0] >= 0, erod_dri.sum(axis=-1), erod_total)
    if dust_soils == 1:
        if clay_nga is None or sand_nga is None:
            raise ValueError("[dust_afwa] dust_soils=1 requires clay_nga and sand_nga")
        available = clay_nga >= 0
        clay = np.where(available, clay_nga, clay)
        sand = np.where(available, sand_nga, sand)
    if dust_veg == 1:
        if vegfra is None:
            raise ValueError("[dust_afwa] dust_veg=1 requires vegfra in percent")
        vegmask = np.where(vegfra >= settings.veg_cutoff, 0.0, 1.0)
    elif dust_veg == 2:
        if lai_vegmask is None:
            raise ValueError("[dust_afwa] dust_veg=2 requires lai_vegmask")
        vegmask = lai_vegmask
    else:
        vegmask = (erod[..., 0] != 0).astype(float)

    return erod_total * vegmask, clay, sand


def _afwa_thresholds(
    volsm, clay, airden, isltyp, dust_smois=0, sf_surface_physics=-1, settings=None
):
    """Shared land-only wet thresholds, including bin seven for lofting."""
    settings = AfwaSettings() if settings is None else settings
    soil = isltyp
    if ((soil < 1) | (soil > len(settings.porosity_table)) | (soil != np.floor(soil))).any():
        raise ValueError("[dust_afwa] land isltyp must be an integer WRF soil type from 1 to 19")
    poros = np.asarray(settings.porosity_table)[soil.astype(int) - 1]
    volsm_t = np.maximum(volsm * settings.smtune, 0.0)
    denom = (1.0 - poros) * (
        settings.soil_dens_sand * (1.0 - clay) + settings.soil_dens_clay * clay
    )
    gravsm = 100.0 * volsm_t / denom
    smois_opt = 0
    if dust_smois == 1 and sf_surface_physics in (3, 7):
        drylimit = (
            settings.dry_ruc_scale
            * (settings.dry_ruc_clay * clay + settings.dry_ruc_offset) ** settings.dry_ruc_power
        )
        smois_opt = 1
    elif dust_smois == 1 and sf_surface_physics == 2:
        drylimit = (
            settings.dry_noah_scale
            * (settings.dry_noah_clay * clay + settings.dry_noah_offset) ** settings.dry_noah_power
        )
        smois_opt = 1
    else:
        drylimit = settings.dry_base_quad * clay**2 + settings.dry_base_linear * clay
    moisture = 100.0 * volsm_t if smois_opt == 1 else gravsm
    u_ts0 = threshold_friction_velocity(settings.reff_salt, settings.den_salt, airden, settings)
    wet = (moisture > drylimit)[:, :, np.newaxis]
    corr = np.sqrt(
        1.0 + settings.moist_scale * np.maximum(moisture - drylimit, 0.0) ** settings.moist_power
    )
    u_ts = np.where(wet, u_ts0 * corr[:, :, np.newaxis], u_ts0)
    return np.maximum(u_ts, 0.0)


def dust_emission_afwa(
    ustar,
    volsm,
    airden,
    znt,
    snowh,
    erod,
    clay,
    sand,
    oro,
    isltyp,
    dust_veg=0,
    vegfra=None,
    lai_vegmask=None,
    dust_smois=0,
    sf_surface_physics=-1,
    dust_dsr=0,
    erod_dri=None,
    dust_soils=0,
    clay_nga=None,
    sand_nga=None,
    settings=None,
):
    """
    Compute AFWA dust emission flux (WRF-Chem source_dust, vectorized).
    Only land cells are evaluated; sea flux is zero (afwa.F L114).

    erod contains the 3 Ginoux bins; Python sums them (afwa.F L161).
    With dust_veg other than 1 or 2 (Registry L3909 default 0), the mask is 0 where Ginoux
    bin 1 is zero, else 1 (afwa.F L185-191). With dust_veg=1, vegfra is
    required in percent: values >=5 are masked (L166-176). With dust_veg=2,
    lai_vegmask is required and used directly (L177-182). Total erodibility
    is multiplied by this internally constructed mask (L196).
    isltyp is required; soil types 15, 16, 18 are always excluded (L238-241).
    Porosity is derived from isltyp (L257); cmb is fixed at 1 (L488).
    dust_dsr=1 selects DRI totals where its first bin is >=0 (L154-162).
    dust_soils=1 selects NGA textures where clay_nga>=0 (L206-217); these
    textures are used in every later calculation. Missing selected datasets
    raise ValueError. The Ginoux mask always uses erod, never erod_dri (L187).
    dust_smois=1 and sf_surface_physics 3 or 7 use the RUC/PX dry limit;
    surface model 2 uses Noah's dry limit (Registry.EM_COMMON L3210,3211,3215).
    Those branches compare 100*volsm with the dry limit; every other case
    uses gravimetric moisture (afwa.F L264-285,558-572).

    Parameters
    ----------
    ustar   : np.ndarray [x, y]   friction velocity u* [m/s]
    volsm   : np.ndarray [x, y]   volumetric soil moisture [m^3/m^3]
    airden  : np.ndarray [x, y]   near-surface air density [kg/m^3]
    znt     : np.ndarray [x, y]   roughness length [m]
    snowh   : np.ndarray [x, y]   snow depth [m]
    erod    : np.ndarray [x, y, 3] Ginoux erodibility fractions [-]
    clay    : np.ndarray [x, y]   clay mass fraction [-]
    sand    : np.ndarray [x, y]   sand mass fraction [-]
    oro     : np.ndarray [x, y]   land-water mask (WRF xland < 1.5 = land)
    isltyp  : np.ndarray [x, y]   soil-type index, excludes 15/16/18 (required)
    dust_veg : int   vegetation selector: 1, 2 or Ginoux fallback (default 0)
    vegfra, lai_vegmask : np.ndarray [x, y]   required for options 1 and 2
    dust_smois, sf_surface_physics : int   moisture and LSM selectors
    dust_dsr, dust_soils : int   regional dataset selectors (default 0)
    erod_dri : np.ndarray [x, y, 3]   required when dust_dsr=1
    clay_nga, sand_nga : np.ndarray [x, y]   required when dust_soils=1
    alpha, gamma, smtune, ustune : float   tuning constants (see header)

    Returns
    -------
    F_p : np.ndarray [x, y, n_dust]   dust emission flux [kg m^-2 s^-1]
    """
    settings = AfwaSettings() if settings is None else settings
    fields = locals().copy()
    shape = np.shape(oro)
    if np.shape(erod) != (*shape, 3):
        raise ValueError("[dust_afwa] erod must have shape [x, y, 3]")
    for name, value in (
        ("ustar", ustar),
        ("volsm", volsm),
        ("airden", airden),
        ("znt", znt),
        ("snowh", snowh),
        ("clay", clay),
        ("sand", sand),
        ("isltyp", isltyp),
    ):
        if np.shape(value) != shape:
            raise ValueError(f"[dust_afwa] {name} shape must match oro")
    chosen = (
        ("erod_dri", erod_dri, (*shape, 3)) if dust_dsr == 1 else None,
        ("clay_nga", clay_nga, shape) if dust_soils == 1 else None,
        ("sand_nga", sand_nga, shape) if dust_soils == 1 else None,
        ("vegfra", vegfra, shape) if dust_veg == 1 else None,
        ("lai_vegmask", lai_vegmask, shape) if dust_veg == 2 else None,
    )
    for item in chosen:
        if item is not None:
            name, value, expected = item
            if value is None or np.shape(value) != expected:
                raise ValueError(
                    f"[dust_afwa] selected option requires {name} with shape {expected}"
                )
    land = oro < 1.5
    if not land.all():
        flux = np.zeros((*oro.shape, NDUST))
        if land.any():
            flux[land] = dust_emission_afwa(**_afwa_land_fields(fields, land))[:, 0, :]
        return flux
    erod_eff, clay, sand = _afwa_surface_fields(
        erod,
        clay,
        sand,
        dust_dsr,
        erod_dri,
        dust_soils,
        clay_nga,
        sand_nga,
        dust_veg,
        vegfra,
        lai_vegmask,
        settings,
    )
    _warn_land_inputs(erod, clay, sand, volsm)

    us = (ustar * settings.ustune)[:, :, np.newaxis]
    u_ts = _afwa_thresholds(volsm, clay, airden, isltyp, dust_smois, sf_surface_physics, settings)

    # --- saltation surface-area weighting -----------------------------------
    ds_rel = saltation_area_fraction(clay, sand, settings)  # [x,y,smx]

    # --- emission mask: land, erodible, not too rough, snow-free ------------
    emit_mask = land & (erod_eff > 0.0) & ~(znt > settings.znt_max) & ~(snowh > settings.snow_max)
    non_erodible = np.isin(isltyp, settings.nonerodible_soil)
    emit_mask = emit_mask & (~non_erodible)

    # --- saltation flux per bin (MB95), only where u* > u_ts ----------------
    can_salt = (us > u_ts) & emit_mask[:, :, np.newaxis]
    # Inactive cells use a safe denominator; emitting cells retain the same arithmetic.
    salt_us = np.where(can_salt, us, 1.0)
    salt = np.where(
        can_salt,
        settings.cmb
        * ds_rel
        * (airden[:, :, np.newaxis] / GRAV)
        * salt_us**3
        * (1.0 + u_ts / salt_us)
        * (1.0 - (u_ts**2) / (salt_us**2)),
        0.0,
    )  # [x,y,smx] kg m^-1 s^-1

    # --- sandblasting efficiency beta (clay-driven, capped) -----------------
    beta = np.minimum(
        10.0 ** (settings.beta_clay * clay + settings.beta_offset), settings.betamax
    )  # [x,y]

    # --- total vertical (emitted) mass flux ---------------------------------
    salt_total = salt.sum(axis=-1)  # [x,y] kg m^-1 s^-1
    emit = salt_total * (erod_eff**settings.gamma) * settings.alpha * beta  # [x,y] kg m^-2 s^-1
    emit = np.where(emit_mask, emit, 0.0)
    emit = np.maximum(emit, 0.0)

    # --- split into dust bins via fixed Kok distribution --------------------
    F_p = emit[:, :, np.newaxis] * np.asarray(settings.distr_dust)[np.newaxis, np.newaxis, :]
    return F_p  # [x, y, n_dust] kg m^-2 s^-1


def afwa_dust_driver(
    dust,
    rho_phy,
    dz8w,
    smois,
    dt,
    u10,
    v10,
    emis_dust,
    tot_edust,
    dustin,
    chem_opt=0,
    settings=None,
    **fields,
):
    """Return the AFWA default-GOCART driver outputs without modifying inputs.

    dust is [x, level, y, 5] in ug/kg; rho_phy is [x, level, y] in
    kg/m3. dz8w is [x, level, y] in m; smois is [x, soil layer, y].
    The driver takes level 1 and soil layer 1. dt is seconds.
    u10/v10, tot_edust are [x, y]; emis_dust and dustin are [x, y, 5].
    emis_dust/tot_edust are accumulated kg/m2. fields are the unchanged
    dust_emission_afwa arguments except airden, taken from rho_phy[:,0,:].
    Returned dust retains all levels; tot_dust/vis_dust are [x, level, y]
    in ug/m3 and m. afwa_dustloft is [x,y] in m/s, -99 on sea or where
    effective erodibility is not positive. flux retains kg/m2/s units.

    chem_opt 2/11 writes the updated lowest tracers to dustin and leaves
    dust unchanged (F L293-294). Other supported options use the final
    GOCART ELSE (F L311-317). CB05_SORG_AQ_KPP and
    CB05_SORG_VBS_AQ_KPP (F L120-126,295-310) are NOT ported: they
    transfer dust to other chemistry packages; do not pass those options.
    Diagnostics follow F L318-367, including every level over sea.
    """
    settings = AfwaSettings() if settings is None else settings
    dz8w = np.asarray(dz8w)[:, 0, :]
    fields = dict(
        fields, airden=rho_phy[:, 0, :], volsm=np.asarray(smois)[:, 0, :], settings=settings
    )
    flux = dust_emission_afwa(**fields)
    land = fields["oro"] < 1.5
    updated = np.array(dust, dtype=float, copy=True)
    emitted = np.array(emis_dust, dtype=float, copy=True)
    total_emitted = np.array(tot_edust, dtype=float, copy=True)
    incoming = np.array(dustin, dtype=float, copy=True)
    loft = np.full(land.shape, -99.0)
    if land.any():
        f = _afwa_land_fields(fields, land)
        timestep = np.broadcast_to(dt, land.shape)[land]
        mass = np.maximum(flux[land] * timestep[:, None], 0.0)
        lowest = updated[:, 0, :, :]
        tc = lowest[land] * 1.0e-9 + mass / dz8w[land, None] / fields["airden"][land, None]
        if chem_opt in (2, 11):
            incoming[land] = tc * 1.0e9
        else:
            lowest[land] = tc * 1.0e9
        emitted[land] = emitted[land] + mass
        total_emitted[land] = total_emitted[land] + mass.sum(axis=-1)
        erod_eff, clay, _ = _afwa_surface_fields(
            f["erod"],
            f["clay"],
            f["sand"],
            **{
                name: f[name]
                for name in (
                    "dust_dsr",
                    "erod_dri",
                    "dust_soils",
                    "clay_nga",
                    "sand_nga",
                    "dust_veg",
                    "vegfra",
                    "lai_vegmask",
                )
                if name in f
            },
            settings=settings,
        )
        threshold = _afwa_thresholds(
            f["volsm"],
            clay,
            f["airden"],
            f["isltyp"],
            f.get("dust_smois", 0),
            f.get("sf_surface_physics", -1),
            f["settings"],
        )[..., 6][:, 0]
        tune = np.broadcast_to(settings.ustune, land.shape)[land]
        us = fields["ustar"][land] * tune
        wind = np.sqrt(u10[land] ** 2 + v10[land] ** 2)
        rough = fields["znt"][land]
        psi = np.zeros_like(wind)
        active = (us != 0.0) & (rough != 0.0)
        # WRF permits nonfinite lofting for singular land inputs (F L322-327).
        with np.errstate(divide="ignore", invalid="ignore"):
            psi[active] = 0.4 * wind[active] / us[active] - np.log(
                settings.loft_height / rough[active]
            )
            potential = tune * wind - threshold * (np.log(settings.loft_height / rough) + psi) / 0.4
        loft[land] = np.where(erod_eff[:, 0] > 0.0, potential, -99.0)
    concentration = updated.sum(axis=-1) * rho_phy
    extinction = (updated * np.asarray(settings.extinction)).sum(axis=-1) * rho_phy
    visibility = np.full(concentration.shape, 999999.0)
    positive = concentration > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        visibility[positive] = np.minimum(
            settings.visibility_scale / extinction[positive], 999999.0
        )
    return dict(
        flux=flux,
        dust=updated,
        afwa_dustloft=loft,
        emis_dust=emitted,
        tot_edust=total_emitted,
        tot_dust=concentration,
        vis_dust=visibility,
        dustin=incoming,
    )


# ---------------------------------------------------------------------------
# Box-model demo
# Run with: python dust_afwa.py
# ---------------------------------------------------------------------------
