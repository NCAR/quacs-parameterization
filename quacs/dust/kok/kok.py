"""
Dust Emission Scheme: Kok et al. (2014) / K14
==============================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Kok, J. F., et al. (2014). An improved dust emission model – Part 1:
    Model description and comparison against measurements.
    Atmos. Chem. Phys., 14, 13023–13041.
    https://doi.org/10.5194/acp-14-13023-2014

    Kok, J. F., et al. (2014). An improved dust emission model – Part 2:
    Evaluation in the Community Earth System Model, with implications for
    the use of dust source functions.
    Atmos. Chem. Phys., 14, 13043–13061.
    https://doi.org/10.5194/acp-14-13043-2014

    Shao, Y. and Lu, H. (2000). A simple expression for wind erosion
    threshold friction velocity. J. Geophys. Res., 105, 22437–22443.
    https://doi.org/10.1029/2000JD900304

    Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization
    of the increase of the aeolian erosion threshold wind friction velocity
    due to soil moisture for arid and semi-arid areas.
    Ann. Geophys., 17, 149–157.
    https://doi.org/10.1007/s00585-999-0149-7

    Fortran reference: CESM2 / CLM5
      DUSTMod.F90       dust emission flux (K14 branch)
      DustEmisLeung2023.F90  updated scheme (drag partition, intermittency)
    https://github.com/ESCOMP/CTSM

QUACS Contract
--------------
Standalone Python module. No MPAS-A, no config files, no external
dependencies beyond NumPy. Import it as ``quacs.dust.kok``.
The unit tests are in ``tests/test_kok.py``. Run them with
``pytest quacs/dust/kok``. The examples are in ``examples/``.

Pre-API note:
    Once the QUACS scripting API is finalized, the expected signature is:
        @scheme(category='dust_emission')
        def dust_emission_k14(state: ModelState, config: Config) -> SurfaceFlux:
    Until then, plain NumPy arrays are used. This script is the scientific
    specification that drives that API design.

Physics summary (Kok et al. 2014, Eq. 7a–7b):

    Standardized threshold friction velocity (proxy for soil aridity):
        u_st = u_ft * sqrt(rho_a / rho_a0)

    Soil erodibility coefficient (exponential in u_st):
        C_d = C_d0 * exp(-C_e * (u_st - u_st0) / u_st0)

    Fragmentation exponent (sensitivity of F_d to u_*):
        kappa = C_kappa * (u_st - u_st0) / u_st0        [capped at KAPPA_MAX]

    Dust emission flux (Kok et al. 2014a; same form in Leung et al. 2023):
        F_d = C_tune * C_d * f_bare * f_clay
              * rho_a * (u_s^2 - u_t^2) / u_st * (u_s / u_t)^kappa
              [only when u_s > u_t]

        Units: [kg m-3] * [m2 s-2] / [m s-1] = [kg m-2 s-1]
        (C_tune, C_d, f_bare, f_clay, and the (u_s/u_t)^kappa term are
        dimensionless.)

    where the "wet" fluid threshold friction velocity u_ft accounts for
    soil moisture following Fécan et al. (1999):
        f_m = sqrt(1 + 1.21 * (100*(w - w_t))^0.68)   [when w > w_t]
        u_ft = u_ft0 * f_m

    and u_ft0 (dry threshold) is computed from Shao & Lu (2000):
        u_ft0 = sqrt(A * (rho_p * g * D_p + gamma / D_p) / rho_a)

    Emission is zeroed over ocean/sea-ice (oro != 1), over snow-covered
    cells (snowdp > SNOW_THRESH), and when VAI exceeds VAI_THR.

Inputs (dynamic, per time step from CLM5):
    ustar  : np.ndarray [x, y]     surface friction velocity [m/s]
    w      : np.ndarray [x, y]     gravimetric soil moisture [-], in [0, 1]
    rho_a  : np.ndarray [x, y]     near-surface air density [kg/m^3]
    vai    : np.ndarray [x, y]     vegetation area index (LAI + SAI) [-]
    snowdp : np.ndarray [x, y]     snow depth [m]

Inputs (static, loaded once at initialization):
    f_clay : np.ndarray [x, y]     clay fraction [-], in [0, 1]
    oro    : np.ndarray [x, y]     land-water mask (1.0=land, 2.0=ocean, etc.)

Configurable:
    C_tune  : float  global tuning constant (dimensionless).
                     Set so that modeled global mean DAOD ≈ 0.03 following
                     Ridley et al. (2016). Kok et al. (2014b) used C_tune=0.05
                     to normalize K14 global emission to match Z03.
                     Default here preserves that published value.
    Dp      : float  median soil particle diameter [m].
                     CLM5 default (Z03): 75e-6 m (optimal lift diameter).
                     Leung et al. (2023) revised value: 127e-6 m (observed
                     median across arid-region soil PSDs).
                     Default here uses the K14/CLM5 value of 75e-6 m.

Output:
    F_d : np.ndarray [x, y]   total surface dust emission flux [kg m^-2 s^-1]
                               (integrated across all emitted size bins;
                               CAM6 distributes into Aitken/accumulation/coarse
                               modes via brittle-fragmentation size distribution,
                               Kok 2011. QUACS runtime converts to tendency:
                               dx/dt = F_d / (rho_air * dz))

Note on the emission threshold u_t:
    In the standard K14 scheme (Kok et al. 2014a,b), the dust emission
    threshold is set to the fluid (static) threshold u_ft. The Leung et al.
    (2023) intermittency extension replaces u_ft with the impact (dynamic)
    threshold u_it = B_it * u_ft0, allowing emissions to continue at lower
    winds once saltation is initiated. This file implements the baseline K14
    physics (u_t = u_ft). The Leung 2023 intermittency scheme is implemented
    separately in quacs/dust/leung.

Note on f_clay capping:
    Following Leung et al. (2023) / CLM5, f_clay in the emission equation
    is capped at FCLAY_MAX = 0.2 (i.e., fclay0 = min(f_clay, 0.2)) to
    prevent unrealistically high fluxes over clay-dominated soils.

Note on VAI threshold:
    K14 (following Mahowald et al. 2010) uses VAI_THR = 0.3, above which
    f_bare = 0 and no dust is emitted. Leung et al. (2023) revised this to
    VAI_THR = 1.0 to allow emissions from semiarid regions. VAI_THR is
    exposed as a module-level constant so callers can override it.
"""

import os

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80616  # gravitational acceleration [m/s^2]
RHO_A0 = 1.225  # standard atmospheric density [kg/m^3]  (ICAO standard)
LAND = 1.0  # oro value identifying land cells

# ---------------------------------------------------------------------------
# K14 emission-equation constants (Kok et al. 2014a, Table 1)
# ---------------------------------------------------------------------------
CD0 = 4.4e-5  # dust emission coefficient prefactor [-]
CE = 2.0  # erodibility decay exponent [-]
U_ST0 = 0.16  # minimum standardized threshold (optimally erodible) [m/s]
C_KAPPA = 2.7  # fragmentation exponent coefficient [-]
KAPPA_MAX = 3.0  # upper cap on kappa (Leung et al. 2023, Sect. 3.4)

# Tuning constant (Kok et al. 2014b used 0.05 to normalize to Z03 global total)
C_TUNE_DEFAULT = 0.05  # dimensionless

# ---------------------------------------------------------------------------
# Shao & Lu (2000) threshold constants
# ---------------------------------------------------------------------------
SL_A = 0.0123  # empirical constant [-]
SL_GAMMA = 1.65e-4  # interparticle force coefficient [kg s^-2]

# Default soil median particle diameter [m]
# K14 / CLM5 default (optimal lift diameter in I&W82): 75 µm
# Leung et al. (2023) observational median: 127 µm
DP_DEFAULT = 75.0e-6  # [m]

# Particle bulk density (quartz, used in Shao & Lu 2000 threshold)
RHOP_DEFAULT = 2650.0  # [kg/m^3]

# ---------------------------------------------------------------------------
# Fécan et al. (1999) soil-moisture correction constants
# ---------------------------------------------------------------------------
# Tuning constant 'a' for threshold gravimetric soil moisture w_t.
# CLM5 / Z03 uses a = 1/f_clay.  K14 (Kok et al. 2014b) uses a = 1.
# Leung et al. (2023) CESM2 evaluation uses a = 2 to compensate for
# CLM5 high-moisture bias.
FECAN_A_DEFAULT = 1.0  # K14 published value

# ---------------------------------------------------------------------------
# Vegetation and land-surface thresholds
# ---------------------------------------------------------------------------
# K14 / Mahowald et al. (2010): dust emission ceases when VAI >= 0.3
# Leung et al. (2023): revised to 1.0 to capture semiarid sources
VAI_THR = 0.3  # [-]; override to 1.0 for Leung 2023 behaviour

# f_clay cap (CLM5 / Leung et al. 2023)
FCLAY_MAX = 0.2  # [-]

# Snow suppression threshold [m]  — emission is zeroed above this depth
SNOW_THRESH = 0.01  # [m]


# ---------------------------------------------------------------------------
# Shao & Lu (2000) dry threshold friction velocity
# ---------------------------------------------------------------------------
def threshold_velocity_dry(Dp=DP_DEFAULT, rhop=RHOP_DEFAULT, rho_a=None):
    """
    Compute the dry fluid threshold friction velocity u_ft0 [m/s] using the
    Shao & Lu (2000) expression (Leung et al. 2023, Eq. 6):

        u_ft0 = sqrt(A * (rho_p * g * D_p  +  gamma / D_p) / rho_a)

    This is computationally simpler than the Iversen & White (1982) scheme
    used in the original CLM5/Z03 and avoids an iterative solution.

    Parameters
    ----------
    Dp   : float or np.ndarray   soil median particle diameter [m].
                                 Scalar or broadcastable to rho_a shape.
    rhop : float                 soil particle bulk density [kg/m^3].
    rho_a : np.ndarray [x, y]   near-surface air density [kg/m^3].
                                 If None, uses RHO_A0 (scalar, standard atm).

    Returns
    -------
    u_ft0 : np.ndarray [x, y]   dry fluid threshold friction velocity [m/s].
            Shape matches rho_a; scalar if rho_a is None.
    """
    if rho_a is None:
        rho_a = RHO_A0
    u_ft0 = np.sqrt(SL_A * (rhop * GRAV * Dp + SL_GAMMA / Dp) / rho_a)
    return u_ft0


# ---------------------------------------------------------------------------
# Fécan et al. (1999) soil-moisture correction factor
# ---------------------------------------------------------------------------
def moisture_factor(w, f_clay, a=FECAN_A_DEFAULT):
    """
    Compute the soil-moisture enhancement factor f_m >= 1 following
    Fécan et al. (1999), as implemented in CLM5 (Leung et al. 2023, Eq. 2):

        w_t = 0.01 * a * (17 * f_clay  +  14 * f_clay^2)   [fraction]
        f_m = sqrt(1 + 1.21 * (100*(w - w_t))^0.68)   when w > w_t
        f_m = 1                                          when w <= w_t

    where w and w_t are both expressed as fractions (kg water / kg soil).

    Parameters
    ----------
    w      : np.ndarray [x, y]  gravimetric soil moisture [kg/kg], in [0, 1].
    f_clay : np.ndarray [x, y]  clay fraction [-], in [0, 1].
    a      : float              tuning constant for w_t (default 1.0 per K14).

    Returns
    -------
    f_m : np.ndarray [x, y]   moisture enhancement factor [-], >= 1.
    """
    w_t = 0.01 * a * (17.0 * f_clay + 14.0 * f_clay**2)  # [fraction]
    excess = np.maximum(100.0 * (w - w_t), 0.0)  # [percent excess]
    f_m = np.where(w > w_t, np.sqrt(1.0 + 1.21 * excess**0.68), 1.0)
    return f_m


# ---------------------------------------------------------------------------
# Standardized threshold and soil erodibility coefficient
# ---------------------------------------------------------------------------
def standardized_threshold(u_ft, rho_a):
    """
    Compute the standardized fluid threshold friction velocity u_st [m/s]
    (Leung et al. 2023, Eq. 5a):

        u_st = u_ft * sqrt(rho_a / rho_a0)

    u_st is a pure function of soil moisture (the sqrt(rho_a) cancels the
    rho_a^{-0.5} dependence in u_ft0 from Shao & Lu 2000). It serves as
    a proxy for soil aridity and drives the erodibility coefficient C_d.

    Parameters
    ----------
    u_ft  : np.ndarray [x, y]   wet fluid threshold [m/s].
    rho_a : np.ndarray [x, y]   near-surface air density [kg/m^3].

    Returns
    -------
    u_st : np.ndarray [x, y]   standardized fluid threshold [m/s].
    """
    return u_ft * np.sqrt(rho_a / RHO_A0)


def erodibility_coefficient(u_st):
    """
    Compute the soil erodibility coefficient C_d (Kok et al. 2014a, Eq. 7b):

        C_d = C_d0 * exp(-C_e * (u_st - u_st0) / u_st0)

    C_d increases exponentially as u_st decreases (drier, more erodible soil).
    This is the key physics that K14 adds over Z03: soil erodibility is not
    a static source function but a dynamic, physically based quantity.

    Parameters
    ----------
    u_st : np.ndarray [x, y]   standardized fluid threshold [m/s].

    Returns
    -------
    C_d : np.ndarray [x, y]   soil erodibility coefficient [-].
    """
    return CD0 * np.exp(-CE * (u_st - U_ST0) / U_ST0)


def fragmentation_exponent(u_st):
    """
    Compute the fragmentation exponent kappa (Kok et al. 2014a, Eq. 7a):

        kappa = C_kappa * (u_st - u_st0) / u_st0

    kappa quantifies how sensitively F_d responds to wind speed above
    threshold. It is capped at KAPPA_MAX = 3.0 following Leung et al. (2023)
    to avoid unrealistically large fluxes at high u_st (low-moisture regions).

    Parameters
    ----------
    u_st : np.ndarray [x, y]   standardized fluid threshold [m/s].

    Returns
    -------
    kappa : np.ndarray [x, y]   fragmentation exponent [-], <= KAPPA_MAX.
    """
    kappa = C_KAPPA * (u_st - U_ST0) / U_ST0
    return np.minimum(kappa, KAPPA_MAX)


# ---------------------------------------------------------------------------
# Bare fraction
# ---------------------------------------------------------------------------
def bare_fraction(vai, vai_thr=VAI_THR):
    """
    Compute the bare-soil fraction available for dust emission.

    Following CLM5 / Mahowald et al. (2010):
        f_bare = max(0,  1 - vai / vai_thr)

    Emission is zero when vai >= vai_thr (default 0.3 for K14;
    override to 1.0 for Leung et al. 2023 behaviour).

    Parameters
    ----------
    vai     : np.ndarray [x, y]   vegetation area index (LAI + SAI) [-].
    vai_thr : float               VAI threshold above which no emission occurs.

    Returns
    -------
    f_bare : np.ndarray [x, y]   bare-soil fraction [-], in [0, 1].
    """
    return np.maximum(0.0, 1.0 - vai / vai_thr)


# ---------------------------------------------------------------------------
# Initialization: load static external datasets once at model start
# ---------------------------------------------------------------------------
def initialize(f_clay_file, oro_file):
    """
    Load and validate static input fields. Called once before the time loop.

    Parameters
    ----------
    f_clay_file : str   Path to clay fraction [x, y], values in [0, 1].
    oro_file    : str   Path to land-water mask [x, y] (1.0 = land).

    Returns
    -------
    f_clay : np.ndarray [x, y]
    oro    : np.ndarray [x, y]
    """
    files = [f_clay_file, oro_file]
    labels = ["f_clay_file", "oro_file"]
    for path, label in zip(files, labels):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_kok2014] {label} not found: {path}")

    f_clay = np.load(f_clay_file)
    oro = np.load(oro_file)

    if f_clay.ndim != 2:
        raise ValueError(f"[dust_kok2014] f_clay must be 2D [x,y], got {f_clay.shape}")
    if not (f_clay.min() >= 0.0 and f_clay.max() <= 1.0):
        raise ValueError(
            f"[dust_kok2014] f_clay must be in [0,1]; got [{f_clay.min():.3f}, {f_clay.max():.3f}]"
        )
    if oro.shape != f_clay.shape:
        raise ValueError(f"[dust_kok2014] oro shape {oro.shape} != f_clay shape {f_clay.shape}")

    return f_clay, oro


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def dust_emission(
    ustar,
    w,
    rho_a,
    vai,
    snowdp,
    f_clay,
    oro,
    Dp=DP_DEFAULT,
    rhop=RHOP_DEFAULT,
    a=FECAN_A_DEFAULT,
    vai_thr=VAI_THR,
    C_tune=C_TUNE_DEFAULT,
):
    """
    Compute total dust emission flux: Kok et al. (2014), Eq. 7a–7b.

    Physics:

        u_ft0 = sqrt(A * (rho_p * g * Dp  +  gamma / Dp) / rho_a)   [S&L00]
        f_m   = Fécan et al. (1999) moisture factor
        u_ft  = u_ft0 * f_m                                          [wet threshold]
        u_t   = u_ft                                                  [K14: u_t = u_ft]
        u_st  = u_ft * sqrt(rho_a / rho_a0)                          [standardized]
        C_d   = C_d0 * exp(-C_e * (u_st - u_st0) / u_st0)           [erodibility]
        kappa = min(C_kappa * (u_st - u_st0)/u_st0,  kappa_max)      [fragmentation]
        f_bare = max(0, 1 - vai / vai_thr)                           [bare fraction]
        fclay0 = min(f_clay, 0.2)                                    [clay cap]

        F_d = C_tune * C_d * f_bare * fclay0
              * rho_a * (ustar^2 - u_t^2) / u_st * (ustar / u_t)^kappa
              [only where ustar > u_t, oro == LAND, snowdp <= SNOW_THRESH]

        Units: [kg m-3] * [m2 s-2] / [m s-1] = [kg m-2 s-1]

    Emission is zero over non-land cells (oro != 1.0), over snow-covered
    cells (snowdp > SNOW_THRESH), when ustar <= u_ft, or when f_bare = 0.

    Parameters
    ----------
    ustar   : np.ndarray [x, y]   surface friction velocity [m/s]
    w       : np.ndarray [x, y]   gravimetric soil moisture [kg/kg], in [0, 1]
    rho_a   : np.ndarray [x, y]   near-surface air density [kg/m^3]
    vai     : np.ndarray [x, y]   vegetation area index (LAI + SAI) [-]
    snowdp  : np.ndarray [x, y]   snow depth [m]
    f_clay  : np.ndarray [x, y]   clay fraction [-], in [0, 1]
    oro     : np.ndarray [x, y]   land-water mask (1.0 = land)
    Dp      : float               soil median particle diameter [m]
    rhop    : float               soil particle density [kg/m^3]
    a       : float               Fécan moisture tuning constant
    vai_thr : float               VAI threshold for bare-fraction calculation
    C_tune  : float               global tuning constant (dimensionless)

    Returns
    -------
    F_d : np.ndarray [x, y]   total dust emission flux [kg m^-2 s^-1]
    """
    # --- dry threshold (Shao & Lu 2000) --------------------------------------
    u_ft0 = threshold_velocity_dry(Dp=Dp, rhop=rhop, rho_a=rho_a)

    # --- moisture correction (Fécan et al. 1999) ----------------------------
    f_m = moisture_factor(w, f_clay, a=a)
    u_ft = u_ft0 * f_m  # wet fluid threshold [m/s]
    u_t = u_ft  # K14: emission threshold == fluid threshold

    # --- standardized threshold and derived quantities ----------------------
    u_st = standardized_threshold(u_ft, rho_a)
    C_d = erodibility_coefficient(u_st)
    kappa = fragmentation_exponent(u_st)

    # --- land surface factors -----------------------------------------------
    f_bare = bare_fraction(vai, vai_thr=vai_thr)
    fclay0 = np.minimum(f_clay, FCLAY_MAX)  # cap at 0.2

    # --- emission mask: land, no snow, wind above threshold, bare soil ------
    emit_mask = (oro == LAND) & (snowdp <= SNOW_THRESH) & (ustar > u_t) & (f_bare > 0.0)

    # --- emission flux (Kok et al. 2014a, Eq. 7a) ---------------------------
    # Avoid division by zero / negative base in (ustar/u_t)^kappa for masked cells.
    u_t_safe = np.where(emit_mask, u_t, 1.0)  # dummy = 1 where masked
    ustar_safe = np.where(emit_mask, ustar, 1.0)

    # rho_a * (u*^2 - u_t^2) / u_st  ->  [kg m-3][m2 s-2]/[m s-1] = [kg m-2 s-1]
    # (No 1/g factor here: that belongs to the Z03/White horizontal saltation
    #  flux ~ rho_a u*^3 / g, not to the K14 vertical dust flux.)
    flux_term = rho_a * (ustar_safe**2 - u_t_safe**2) / u_st
    shape_term = (ustar_safe / u_t_safe) ** kappa  # [-]

    F_d = np.where(
        emit_mask,
        C_tune * C_d * f_bare * fclay0 * flux_term * shape_term,
        0.0,
    )

    return F_d  # [kg m^-2 s^-1]
