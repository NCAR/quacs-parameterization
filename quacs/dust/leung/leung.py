"""
Dust Emission Scheme: Leung et al. (2023/2024) — CESM2/CLM5 implementation
============================================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Leung, D. M., Kok, J. F., Li, L., Okin, G. S., Prigent, C., Klose, M.,
    Pérez García-Pando, C., Menut, L., Mahowald, N. M., Lawrence, D. M.,
    and Chamecki, M. (2023). A new process-based and scale-aware desert dust
    emission scheme for global climate models – Part I: Description and
    evaluation against inverse modeling emissions.
    Atmos. Chem. Phys., 23, 6487–6523.
    https://doi.org/10.5194/acp-23-6487-2023

    Leung, D. M., Kok, J. F., Li, L., Mahowald, N. M., Lawrence, D. M.,
    Tilmes, S., Kluzek, E., Klose, M., and Pérez García-Pando, C. (2024).
    A new process-based and scale-aware desert dust emission scheme for
    global climate models – Part II: Evaluation in the Community Earth
    System Model version 2 (CESM2).
    Atmos. Chem. Phys., 24, 2287–2318.
    https://doi.org/10.5194/acp-24-2287-2024

    Kok, J. F., Mahowald, N. M., Fratini, G., et al. (2014b). An improved
    dust emission model – Part 2: Evaluation in the Community Earth System
    Model. Atmos. Chem. Phys., 14, 13043–13061.
    https://doi.org/10.5194/acp-14-13043-2014

    Shao, Y. and Lu, H. (2000). A simple expression for wind erosion
    threshold friction velocity. J. Geophys. Res.-Atmos., 105, 22437–22443.
    https://doi.org/10.1029/2000JD900304

    Comola, F., Kok, J. F., Chamecki, M., and Martin, R. L. (2019). The
    Intermittency of Wind-Driven Sand Transport. Geophys. Res. Lett., 46,
    13430–13440. https://doi.org/10.1029/2019GL085739

    Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric
    dust cycle: 1. Design of a soil-derived dust emission scheme.
    J. Geophys. Res., 100, 16415–16430.
    https://doi.org/10.1029/95JD00690

    Okin, G. S. (2008). A new model of wind erosion in the presence of
    vegetation. J. Geophys. Res.-Earth, 113.
    https://doi.org/10.1029/2007JF000758

    Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization
    of the increase of the aeolian erosion threshold wind friction velocity
    due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17,
    149–157. https://doi.org/10.1007/s00585-999-0149-7

    Fortran reference: CLM5 / CESM2 / CTSM
      src/biogeochem/DustEmisLeung2023.F90  (Leung et al. 2023/2024 scheme)
    https://github.com/ESCOMP/CTSM
    https://doi.org/10.5281/zenodo.10621844

QUACS Contract
--------------
Standalone Python module. No MPAS-A, no config files, no external
dependencies beyond NumPy and SciPy. Import it as ``quacs.dust.leung``.
The unit tests are in ``tests/test_leung.py``. Run them with
``pytest quacs/dust/leung``. The examples are in ``examples/``.

Physics summary (Leung et al. 2023/2024):

  Step 1 — Mobilizable bare-ground fraction (fbare):
      fbare = max(0, 1 - VAI/VAIthr) * (1 - frac_sno)
      VAI    = tlai + tsai  (vegetation area index)
      VAIthr = 1.0          (extended from 0.3 in Z03/K14 to capture
                             semiarid marginal sources; Leung et al. 2023)
      Emission occurs only over soil/crop land types with fbare > 0.

  Step 2 — Soil moisture threshold factor (Fécan et al. 1999, Eq. 2):
      gwc_sfc  = h2osoi_vol * rho_water / bd   [kg kg-1]
      bd       = (1 - watsat) * rho_quartz
      gwc_thr  = dust_moist_fact * (0.17 + 0.14*fclay)          [kg kg-1]
                 dust_moist_fact = 1.0
      fm       = sqrt(1 + 1.21*(100*(gwc - gwc_thr))^0.68)  if gwc > gwc_thr
               = 1.0                                         otherwise

      CTSM note (SoilStateInitTimeConstMod.F90): for BOTH the Zender_2003
      and Leung_2023 dust options, CTSM sets
          gwc_thr = dust_moist_fact * ThresholdSoilMoistZender2003(clay)
                  = 1.0 * (0.17 + 0.14 * clay[%] * 0.01)
      i.e. the Zender (2003a) form (Fécan w' x a with a = 1/fclay), not the
      Fécan/K14 form a*(0.17*fclay + 0.14*fclay^2).  An earlier version of
      this file used a = 2 with the Fécan form (0.0368 kg/kg at 10 % clay);
      CTSM gives 0.184 kg/kg at 10 % clay.

  Step 3 — Dry fluid threshold friction velocity (Shao & Lu 2000, Eq. 6):
      u*ft0 = sqrt(A*(rho_p*g*Dp + gamma/Dp) / rho_air)
      Dp    = 127e-6 m  (global median soil diameter; Leung et al. 2023)
      A     = 0.0123, gamma = 1.65e-4 kg s-2
      Units: [kg m-3][m s-2][m] = [kg s-2 m-1]/[kg m-3] -> sqrt -> [m s-1]

  Step 4 — Wet fluid threshold and standardised threshold:
      u*ft  = u*ft0 * fm                        [wet fluid threshold]
      u*st  = u*ft * sqrt(rho_air / rho_air0)   [standardised threshold]

  Step 5 — Impact (dynamic) threshold (Comola et al. 2019, Eq. 11a):
      u*it  = Bit * u*ft0_dry                   Bit = 0.82
      (no soil-moisture factor, as in CTSM)

  Step 6 — Drag partition (Leung et al. 2023, Eqs. 8-10):
      Rock component (Marticorena & Bergametti 1995):
          feff_r = 1 - ln(z0a/z0s) / ln(b1*(X/z0s)^b2)
          z0s = 2*Dp/30, b1=0.7, b2=0.8, X=10 m
      Vegetation component (Okin 2008 / Pierre et al. 2014, Eq. 9):
          K      = 2*(1/fv - 1),  fv = VAI/VAIthr
          feff_v = (K + f0*c) / (K + c),   f0=0.32, c=4.8
      Hybrid drag partition (Eq. 10b):
          Feff^3 = Ar * feff_r^3 + Av * feff_v^3
      Effective soil-surface friction velocity:
          u*s = u* * Feff

  Step 7 — Intermittency factor (Comola et al. 2019; CTSM DustEmisLeung2023):
      All quantities are first converted to *wind speeds* at the saltation
      height z_sal = 0.1 m with the log law (z0 = z0a_glob = 1e-4 m):
          U      = u*s   / kappa_vk * ln(z_sal / z0)      [m s-1]
          U_ft   = u*ft  / kappa_vk * ln(z_sal / z0)      [m s-1]
          U_it   = u*it  / kappa_vk * ln(z_sal / z0)      [m s-1]
      Standard deviation of the wind speed at z_sal (Panofsky et al. 1977):
          sigma_U = u*s * (12 - 0.5 * z_i / L)^(1/3)       [m s-1]
          (z_i = 1000 m; the bracket is floored at 0.001)
      Probabilities and threshold-crossing rate:
          P_ft  = Phi((U_ft - U) / sigma_U)
          P_it  = Phi((U_it - U) / sigma_U)
          alpha = [exp((U_ft^2 - U_it^2 - 2 U (U_ft - U_it)) / (2 sigma_U^2)) + 1]^-1
      Intermittency factor [-]:
          eta = 1 - P_ft + alpha * (P_ft - P_it)

      Mean and sigma are both wind speeds at the same height, so the
      ratio sigma_U / U is dimensionally and physically consistent.

  Step 8 — Soil erodibility coefficient (Kok et al. 2014b, Eq. 5b):
      Cd = Cd0 * exp(-Ce * (u*st - u*st0) / u*st0)
      Cd0=4.4e-5, Ce=2.0, u*st0=0.16 m/s

  Step 9 — Fragmentation exponent (Kok et al. 2014b):
      kappa_raw = Ck * (u*st - u*st0) / u*st0    Ck=2.7
      kappa     = min(kappa_raw, KAPPA_MAX)

  Step 10 — Dust emission flux (Leung et al. 2023; CTSM DustEmisLeung2023):
      if u*s > u*it:
          Fd = eta * Ctune * Cd * fbare * fclay' * liqfrac
               * rho_air * (u*s^2 - u*it^2) / u*it
               * (u*s / u*it)^kappa
      where fclay' = 0.1 + 0.5 * min(fclay, 0.20)        [-], in [0.1, 0.2]
      (CTSM MassFracClayLeung2023 = 0.1 + MassFracClay(clay) * 0.1/0.20,
       MassFracClay = min(clay[%]*0.01, 0.20).  This damps the clay
       dependence relative to K14's min(fclay, 0.2): a clay-free soil still
       has fclay' = 0.1.)

      Units: [kg m-3] * [m2 s-2] / [m s-1] = [kg m-2 s-1]
      (eta, Ctune, Cd, fbare, fclay', liqfrac and (u*s/u*it)^kappa are
      dimensionless.  There is NO 1/g factor: 1/g belongs to the Z03/White
      horizontal saltation flux ~ rho u*^3 / g, not to this vertical flux.)

  Step 11 — Apply upscaling correction map (Leung et al. 2024, Eq. 12):
      Fd_corrected = Fd * Kc_map   (time-invariant spatial scaling)

  Step 12 — Partition into size bins via overlap of lognormal source modes
      with the four transport bins [0.1-1, 1-2.5, 2.5-5, 5-10 µm]:
      F_p[n] = sum_m( ovr_src_snk_mss[m,n] * Fd_corrected )

Inputs (dynamic, per time step from CLM5/MPAS-A):
    fv          : np.ndarray [x, y]    friction velocity [m/s]
    forc_rho    : np.ndarray [x, y]    near-surface air density [kg/m³]
    h2osoi_vol  : np.ndarray [x, y]    volumetric soil water content [m³/m³]
    h2osoi_liq  : np.ndarray [x, y]    liquid soil water [kg/m²]
    h2osoi_ice  : np.ndarray [x, y]    frozen soil water [kg/m²]
    watsat      : np.ndarray [x, y]    saturated volumetric soil water [m³/m³]
    tlai        : np.ndarray [x, y]    one-sided leaf area index [-]
    tsai        : np.ndarray [x, y]    one-sided stem area index [-]
    frac_sno    : np.ndarray [x, y]    snow cover fraction [-]
    obu         : np.ndarray [x, y]    Obukhov length [m] (for intermittency)

Inputs (static, loaded once at initialization):
    fclay           : np.ndarray [x, y]   clay mass fraction [-] (raw; see Steps 2, 10)
    z0a             : np.ndarray [x, y]   aeolian roughness length [m] (Prigent 2005)
    frac_rock       : np.ndarray [x, y]   fractional rock area Ar [-]
    frac_veg        : np.ndarray [x, y]   fractional vegetation area Av [-]
    Kc_map          : np.ndarray [x, y]   upscaling correction factor (Leung 2024)
    is_soil         : np.ndarray [x, y]   boolean mask: soil or crop land type

Precomputed constants (computed once at initialization):
    ovr_src_snk_mss : np.ndarray [3, 4]  source-to-sink bin overlap [-]
    u_ft0_dry_ref   : float               reference dry threshold (at rho_air0)

Output:
    F_p : np.ndarray [x, y, ndst]   surface dust emission flux [kg m-2 s-1]
                                    (ndst = 4 bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm)
                                    QUACS runtime converts to tendency:
                                    dx/dt = F_p / (rho_air * dz)   [kg kg-1 s-1]
                                    or volumetric source:
                                    dC/dt = F_p / dz               [kg m-3 s-1]

Key differences from the DEAD/Z03 scheme:
    1. Base emission equation: K14 (Kok et al. 2014b) replaces White (1979)
       saltation → more physical sensitivity to u* and moisture.
    2. Threshold scheme: Shao & Lu (2000) replaces Iversen & White (1982);
       uses global median Dp = 127 µm instead of optimal-saltation 75 µm.
    3. Threshold used in emission equation: u*it (impact/dynamic) instead of
       u*ft (fluid/static), capturing intermittent emission in marginal regions.
    4. Drag partition: explicit rock + vegetation reduction of u* at the soil
       surface via satellite-derived z0a and VAI-based vegetation scheme.
    5. Emission intermittency: Comola et al. (2019) Gaussian turbulence model
       produces eta ∈ [0,1] within each model time step.
    6. Extended VAI threshold: VAIthr = 1.0 (vs 0.3 in Z03/K14) to capture
       semiarid marginal-source emissions.
    7. Soil moisture threshold: CTSM uses the Zender (2003a) threshold
       gwc_thr = 0.17 + 0.14*fclay (dust_moist_fact = 1) for this scheme,
       and the clay term fclay' = 0.1 + 0.5*min(fclay, 0.2).
    8. Spatial upscaling correction: time-invariant Kc_map adjusts coarse-
       resolution emissions toward the variability of higher-resolution runs.
    9. No preferential source mask S (present in Z03): soil erodibility is
       represented physically through Cd(u*st) instead.
"""

import numpy as np
from scipy.special import erf
from scipy.stats import norm

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80616  # [m s-2]  gravitational acceleration
RHO_WATER = 1000.0  # [kg m-3] density of liquid water
RHO_QUARTZ = 2700.0  # [kg m-3] density of dry soil mineral (quartz)
RHO_AIR0 = 1.225  # [kg m-3] reference air density (sea-level, 15°C)
VON_KARMAN = 0.4  # [-]      von Kármán constant

# ---------------------------------------------------------------------------
# Soil / particle parameters (Leung et al. 2023, Sect. 3.2)
# ---------------------------------------------------------------------------
DP_MED = 127.0e-6  # [m]   global median soil particle diameter
RHO_P = 2650.0  # [kg m-3] soil particle density (quartz)

# Shao & Lu (2000) threshold constants (Eq. 6)
SL_A = 0.0123  # [-]   empirical constant
SL_GAMMA = 1.65e-4  # [kg s-2] interparticle force constant

# ---------------------------------------------------------------------------
# Soil moisture threshold (Fécan et al. 1999 with Zender 2003a tuning; CTSM)
#   gwc_thr = DUST_MOIST_FACT * (0.17 + 0.14 * fclay)        [kg kg-1]
#   CTSM SoilStateInitTimeConstMod.F90 uses this for the Leung_2023 option
#   too, with dust_moist_fact = 1.0.
# ---------------------------------------------------------------------------
DUST_MOIST_FACT = 1.0  # [-]   CTSM dust_moist_fact
FM_A = DUST_MOIST_FACT  # backward-compatible alias (old name)

# ---------------------------------------------------------------------------
# Clay term in the emission equation (CTSM MassFracClayLeung2023)
#   fclay' = FCLAY_LEUNG_MIN + (FCLAY_LEUNG_MIN / FCLAY_MAX) * min(fclay, FCLAY_MAX)
#          = 0.1 + 0.5 * min(fclay, 0.2)                     [-], in [0.1, 0.2]
# ---------------------------------------------------------------------------
FCLAY_MAX = 0.20  # [-]   clay cap (CTSM MassFracClay)
FCLAY_LEUNG_MIN = 0.10  # [-]   offset/scale in MassFracClayLeung2023

# ---------------------------------------------------------------------------
# Impact (dynamic) threshold factor (Comola et al. 2019, Eq. 11a)
# ---------------------------------------------------------------------------
BIT = 0.82  # [-]   u*it = Bit * u*ft0_dry

# ---------------------------------------------------------------------------
# Drag partition constants (Leung et al. 2023, Eq. 8)
# ---------------------------------------------------------------------------
Z0S_FACTOR = 2.0 / 30.0  # [-]   z0s = Z0S_FACTOR * Dp  (Sherman 1992)
DP_B1 = 0.7  # [-]   Marticorena & Bergametti 1995
DP_B2 = 0.8  # [-]   Marticorena & Bergametti 1995
DP_X_ROCK = 10.0  # [m]   IBL length scale for rocks

# Vegetation drag partition (Okin 2008 / Pierre et al. 2014, Eq. 9)
FEFF_V_F0 = 0.32  # [-]   minimum drag partition factor for veg
FEFF_V_C = 4.8  # [-]   constant c in vegetation scheme

# ---------------------------------------------------------------------------
# Vegetation area threshold (Leung et al. 2023; extended from 0.3 to 1.0)
# ---------------------------------------------------------------------------
VAI_MBL_THR = 1.0  # [m² m-2] VAI threshold for dust mobilisation

# ---------------------------------------------------------------------------
# Intermittency (Comola et al. 2019; CTSM DustEmisLeung2023.F90)
# ---------------------------------------------------------------------------
Z_SAL = 0.1  # [m]   saltation height (CTSM hgt_sal)
Z0A_GLOB = 1.0e-4  # [m]   roughness for log-law wind at Z_SAL (CTSM z0a_glob)
ZII = 1000.0  # [m]   convective boundary-layer depth (CTSM zii)
STBLTY_FLOOR = 0.001  # [-]   floor on (12 - 0.5 zii/L) (CTSM)

# ---------------------------------------------------------------------------
# Soil erodibility (Kok et al. 2014b, Eq. 5b)
# ---------------------------------------------------------------------------
CD0 = 4.4e-5  # [-]   reference erodibility coefficient
CE = 2.0  # [-]   erodibility exponential sensitivity
U_ST0 = 0.16  # [m/s] reference standardised threshold

# ---------------------------------------------------------------------------
# Fragmentation exponent (Kok et al. 2014b / Leung et al. 2023, capped)
# ---------------------------------------------------------------------------
CK = 2.7  # [-]   kappa prefactor
KAPPA_MAX = 3.0  # [-]   cap on fragmentation exponent
#       NOTE: CTSM DustEmisLeung2023.F90 uses
#       frag_expt_thr = 2.5; inactive in the demo
#       (kappa ≈ 0.93) but matters for high u*st.

# ---------------------------------------------------------------------------
# Emission flux proportionality (Kok et al. 2014b)
# ---------------------------------------------------------------------------
CTUNE = 0.05  # [-]   global tuning constant

# ---------------------------------------------------------------------------
# Transport size bin boundaries [m]  (ndst = 4)
# Bins: [0.1-1], [1-2.5], [2.5-5], [5-10]  µm
# ---------------------------------------------------------------------------
NDST = 4
DMT_GRD = np.array([0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6])  # [m] edges

# ---------------------------------------------------------------------------
# Source mode parameters (Kok 2011 / BSM96, same lognormal modes as DEAD)
# ---------------------------------------------------------------------------
DST_SRC_NBR = 3
DMT_VMA_SRC = np.array([0.832e-6, 4.82e-6, 19.38e-6])  # [m]  mass median diam
GSD_ANL_SRC = np.array([2.10, 1.90, 1.60])  # [-]  geometric std dev
MSS_FRC_SRC = np.array([0.036, 0.957, 0.007])  # [-]  mass fractions


# ---------------------------------------------------------------------------
# Precomputed constants (computed once at module import)
# ---------------------------------------------------------------------------
def precompute_constants():
    """
    Compute the quantities that are fixed for the lifetime of a simulation:

      1. u_ft0_dry_ref — the "dry" fluid threshold friction velocity at the
                         reference air density RHO_AIR0, using the Shao & Lu
                         (2000) parameterisation at Dp = DP_MED = 127 µm.
                         (Leung et al. 2023, Eq. 6)

      2. ovr_src_snk_mss[m, n] — mass overlap of each source lognormal
                         mode m with each transport sink bin n, used to
                         partition the total vertical dust flux into the four
                         CESM bins.  (identical to DEAD; DUSTMod L777-L785)

    Returns
    -------
    u_ft0_dry_ref   : float             [m/s]  dry threshold at RHO_AIR0
    ovr_src_snk_mss : np.ndarray [3,4]  [-]    source-to-sink mass overlaps
    """
    # --- u_ft0_dry_ref: Shao & Lu (2000) Eq. 6 at rho_air = RHO_AIR0 ------
    numerator = SL_A * (RHO_P * GRAV * DP_MED + SL_GAMMA / DP_MED)
    u_ft0_dry_ref = np.sqrt(numerator / RHO_AIR0)  # [m/s]

    # --- ovr_src_snk_mss: lognormal overlap (same as DEAD precompute) ------
    ovr_src_snk_mss = np.zeros((DST_SRC_NBR, NDST))
    for m in range(DST_SRC_NBR):
        sqrt2lngsd = np.sqrt(2.0) * np.log(GSD_ANL_SRC[m])
        for n in range(NDST):
            lndmax = np.log(DMT_GRD[n + 1] / DMT_VMA_SRC[m])
            lndmin = np.log(DMT_GRD[n] / DMT_VMA_SRC[m])
            ovr_src_snk_frc = 0.5 * (erf(lndmax / sqrt2lngsd) - erf(lndmin / sqrt2lngsd))
            ovr_src_snk_mss[m, n] = ovr_src_snk_frc * MSS_FRC_SRC[m]

    return u_ft0_dry_ref, ovr_src_snk_mss


# Compute module-level constants once at import
U_FT0_DRY_REF, OVR_SRC_SNK_MSS = precompute_constants()


# ---------------------------------------------------------------------------
# Helper: Shao & Lu (2000) dry threshold at arbitrary air density
# ---------------------------------------------------------------------------
def sl2000_u_ft0_dry(forc_rho):
    """
    Dry fluid threshold friction velocity from Shao & Lu (2000) Eq. 6.

    u*ft0 = sqrt( A * (rho_p * g * Dp + gamma / Dp) / rho_air )

    Parameters
    ----------
    forc_rho : np.ndarray [x, y]   near-surface air density [kg/m³]

    Returns
    -------
    u_ft0_dry : np.ndarray [x, y]  [m/s]
    """
    numerator = SL_A * (RHO_P * GRAV * DP_MED + SL_GAMMA / DP_MED)
    return np.sqrt(numerator / forc_rho)


# ---------------------------------------------------------------------------
# Helper: soil moisture threshold factor (Fécan et al. 1999 / Leung 2024)
# ---------------------------------------------------------------------------
def moisture_threshold_factor(h2osoi_vol, watsat, fclay):
    """
    Compute fm, the soil-moisture enhancement factor for u*ft (Eq. 2).

        gwc_sfc = h2osoi_vol * rho_water / bd                    [kg/kg]
        bd      = (1 - watsat) * rho_quartz                       [kg/m³]
        gwc_thr = DUST_MOIST_FACT * (0.17 + 0.14*fclay)           [kg/kg]
        fm      = sqrt(1 + 1.21*(100*(gwc-gwc_thr))^0.68)  if gwc > gwc_thr
                = 1.0                                        otherwise

    gwc_thr follows CTSM (SoilStateInitTimeConstMod.F90), which for the
    Leung_2023 option sets
        gwc_thr = dust_moist_fact * ThresholdSoilMoistZender2003(clay)
                = 1.0 * (0.17 + 0.14 * clay[%] * 0.01)
    This is Fécan's w' = 0.17*fclay + 0.14*fclay^2 multiplied by the
    Zender (2003a) tuning a = 1/fclay.  (Fécan's original form with a = 1
    is what K14 uses; an earlier version of this file used it with a = 2.)

    The factor 100 in fm converts the excess (kg/kg) to percent, as in
    Fécan et al. (1999).

    Parameters
    ----------
    h2osoi_vol : np.ndarray [x, y]  volumetric soil moisture [m³/m³]
    watsat     : np.ndarray [x, y]  saturated volumetric soil water [m³/m³]
    fclay      : np.ndarray [x, y]  clay mass fraction [-]

    Returns
    -------
    fm      : np.ndarray [x, y]  moisture enhancement factor [-], >= 1
    gwc_sfc : np.ndarray [x, y]  surface gravimetric water content [kg/kg]
    gwc_thr : np.ndarray [x, y]  threshold gravimetric water content [kg/kg]
    """
    bd = (1.0 - watsat) * RHO_QUARTZ
    gwc_sfc = h2osoi_vol * RHO_WATER / bd
    gwc_thr = DUST_MOIST_FACT * (0.17 + 0.14 * fclay)  # CTSM Zender2003 threshold
    excess = np.maximum(gwc_sfc - gwc_thr, 0.0)
    fm = np.where(
        gwc_sfc > gwc_thr,
        np.sqrt(1.0 + 1.21 * (100.0 * excess) ** 0.68),
        1.0,
    )
    return fm, gwc_sfc, gwc_thr


# ---------------------------------------------------------------------------
# Helper: clay term used in the emission equation (CTSM)
# ---------------------------------------------------------------------------
def clay_mass_fraction_leung(fclay):
    """
    Clay term fclay' in the Leung emission equation, following CTSM
    MassFracClayLeung2023:

        MassFracClay           = min(clay[%] * 0.01, 0.20)
        MassFracClayLeung2023  = 0.1 + MassFracClay * 0.1 / 0.20

    With fclay as a mass fraction:

        fclay' = 0.1 + 0.5 * min(fclay, 0.20)        [-], in [0.1, 0.2]

    Compared with K14's min(fclay, 0.2), this halves the clay sensitivity
    and gives clay-free soils a non-zero value (0.1).

    Parameters
    ----------
    fclay : np.ndarray [x, y]  clay mass fraction [-]

    Returns
    -------
    fclay_eff : np.ndarray [x, y]  [-]
    """
    return FCLAY_LEUNG_MIN + (FCLAY_LEUNG_MIN / FCLAY_MAX) * np.minimum(fclay, FCLAY_MAX)


# ---------------------------------------------------------------------------
# Helper: drag partition factor (Leung et al. 2023, Eqs. 8-10)
# ---------------------------------------------------------------------------
def drag_partition(z0a, vai, frac_rock, frac_veg):
    """
    Compute the hybrid drag partition factor Feff combining rock (M&B 1995)
    and vegetation (Okin 2008 / Pierre et al. 2014) contributions.

    Rock component (Eq. 8):
        z0s    = Z0S_FACTOR * Dp
        feff_r = 1 - ln(z0a/z0s) / ln(b1 * (X/z0s)^b2)
        clipped to [0, 1]

    Vegetation component (Eq. 9):
        fv     = VAI / VAIthr
        K      = 2 * (1/fv - 1)     (normalized mean gap length)
        feff_v = (K + f0*c) / (K + c)
        Only applied where VAI <= VAIthr; feff_v = 1 where no vegetation.

    Hybrid (Eq. 10b):
        Feff^3 = Ar * feff_r^3 + Av * feff_v^3

    Parameters
    ----------
    z0a      : np.ndarray [x, y]  aeolian roughness length [m] (Prigent 2005)
    vai      : np.ndarray [x, y]  vegetation area index (= tlai + tsai) [-]
    frac_rock: np.ndarray [x, y]  fractional rock area Ar [-]
    frac_veg : np.ndarray [x, y]  fractional vegetation area Av [-]

    Returns
    -------
    Feff   : np.ndarray [x, y]  hybrid drag partition factor [-], in (0, 1]
    feff_r : np.ndarray [x, y]  rock component [-]
    feff_v : np.ndarray [x, y]  vegetation component [-]
    """
    # Smooth roughness length for soil grains (Sherman 1992)
    z0s = Z0S_FACTOR * DP_MED  # scalar [m]

    # ---- Rock drag partition ------------------------------------------------
    ln_ratio = np.log(np.maximum(z0a, 1.0e-10) / z0s)
    ln_denom = np.log(DP_B1 * (DP_X_ROCK / z0s) ** DP_B2)
    feff_r_raw = 1.0 - ln_ratio / ln_denom
    feff_r = np.clip(feff_r_raw, 0.0, 1.0)

    # ---- Vegetation drag partition ------------------------------------------
    fv_safe = np.maximum(vai / VAI_MBL_THR, 1.0e-6)
    K = 2.0 * (1.0 / fv_safe - 1.0)
    feff_v_raw = (K + FEFF_V_F0 * FEFF_V_C) / (K + FEFF_V_C)
    # Apply only where VAI > 0 and VAI <= VAIthr; elsewhere feff_v = 1
    feff_v = np.where(
        (vai > 0.0) & (vai <= VAI_MBL_THR),
        np.clip(feff_v_raw, FEFF_V_F0, 1.0),
        1.0,
    )

    # ---- Hybrid (Eq. 10b) ---------------------------------------------------
    # When no rock or vegetation fractions are assigned (pure bare erodible soil)
    # the drag partition has no roughness elements → Feff = 1.
    total_frac = frac_rock + frac_veg
    Feff3_raw = frac_rock * feff_r**3 + frac_veg * feff_v**3
    Feff3_raw = np.maximum(Feff3_raw, 0.0)
    Feff = np.where(
        total_frac > 0.0,
        np.where(Feff3_raw > 0.0, Feff3_raw ** (1.0 / 3.0), 0.0),
        1.0,
    )

    return Feff, feff_r, feff_v


# ---------------------------------------------------------------------------
# Helper: log-law conversion friction velocity -> wind speed at Z_SAL
# ---------------------------------------------------------------------------
def _log_law_factor():
    """ln(Z_SAL / Z0A_GLOB) / kappa_vk  [-]; U(Z_SAL) = u* x this factor."""
    return np.log(Z_SAL / Z0A_GLOB) / VON_KARMAN


# ---------------------------------------------------------------------------
# Helper: intermittency factor eta (Comola et al. 2019; CTSM)
# ---------------------------------------------------------------------------
def intermittency_factor(u_s, sigma_U, u_it, u_ft):
    """
    Fraction of time within a model time step during which dust emission
    is active, following Comola et al. (2019) as implemented in CTSM
    (DustEmisLeung2023.F90).

    The mean friction velocity u_s and the thresholds u_it, u_ft are
    converted to wind speeds at the saltation height Z_SAL with the log law
    (roughness Z0A_GLOB), so that they can be compared with sigma_U, which
    is the standard deviation of the *wind speed* at Z_SAL:

        U    = u_s  / kappa_vk * ln(Z_SAL / Z0A_GLOB)          [m/s]
        U_ft = u_ft / kappa_vk * ln(Z_SAL / Z0A_GLOB)          [m/s]
        U_it = u_it / kappa_vk * ln(Z_SAL / Z0A_GLOB)          [m/s]

        P_ft  = Phi((U_ft - U) / sigma_U)    (prob. U_tilde < fluid threshold)
        P_it  = Phi((U_it - U) / sigma_U)    (prob. U_tilde < impact threshold)
        alpha = [exp((U_ft^2 - U_it^2 - 2 U (U_ft - U_it)) / (2 sigma_U^2)) + 1]^-1

        eta = 1 - P_ft + alpha * (P_ft - P_it)

    Limits:
      eta -> 1  when U - sigma_U >> U_ft   (continuous emission)
      eta -> 0  when U + sigma_U << U_it   (no emission)

    Parameters
    ----------
    u_s      : np.ndarray [x, y]  mean soil-surface friction velocity [m/s]
    sigma_U  : np.ndarray [x, y]  std dev of wind speed at Z_SAL [m/s]
                                  (from turbulent_sigma_us)
    u_it     : np.ndarray [x, y]  impact (dynamic) threshold friction velocity [m/s]
    u_ft     : np.ndarray [x, y]  fluid (static, wet) threshold friction velocity [m/s]

    Returns
    -------
    eta : np.ndarray [x, y]  intermittency factor in [0, 1]
    """
    lf = _log_law_factor()
    U = u_s * lf
    U_ft = u_ft * lf
    U_it = u_it * lf
    sig = np.maximum(sigma_U, 1.0e-6)

    P_ft = norm.cdf((U_ft - U) / sig)
    P_it = norm.cdf((U_it - U) / sig)

    arg = (U_ft**2 - U_it**2 - 2.0 * U * (U_ft - U_it)) / (2.0 * sig**2)
    alpha = 1.0 / (np.exp(np.clip(arg, -700.0, 700.0)) + 1.0)

    eta = 1.0 - P_ft + alpha * (P_ft - P_it)
    return np.clip(eta, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Helper: turbulent wind std at saltation height (Panofsky et al. 1977)
# ---------------------------------------------------------------------------
def turbulent_sigma_us(u_s, obu):
    """
    Standard deviation of the instantaneous *wind speed* at the saltation
    height Z_SAL = 0.1 m, following Panofsky et al. (1977) as used in
    Comola et al. (2019) and CTSM (DustEmisLeung2023.F90):

        sigma_U = u_s * (12 - 0.5 * z_i / L)^(1/3)          [m/s]

    with z_i = ZII = 1000 m.  The bracket is floored at STBLTY_FLOOR (CTSM),
    which is reached for strongly stable conditions (0 < L < ~42 m).
    Under neutral conditions (|L| -> inf), sigma_U / u_s = 12^(1/3) ≈ 2.29.

    Note: sigma_U is a wind-speed std dev, NOT a friction-velocity std dev.
    It must be compared with wind speeds at Z_SAL (see intermittency_factor).

    Parameters
    ----------
    u_s : np.ndarray [x, y]  soil-surface friction velocity [m/s]
    obu : np.ndarray [x, y]  Obukhov length L [m]; negative = unstable

    Returns
    -------
    sigma_U : np.ndarray [x, y]  [m/s]
    """
    # Avoid division by zero at L = 0 while keeping the sign of L
    obu_safe = np.where(np.abs(obu) < 1.0e-3, np.where(obu < 0.0, -1.0e-3, 1.0e-3), obu)
    bracket = np.maximum(12.0 - 0.5 * ZII / obu_safe, STBLTY_FLOOR)
    return u_s * bracket ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def dust_emission(
    fv,
    forc_rho,
    h2osoi_vol,
    h2osoi_liq,
    h2osoi_ice,
    watsat,
    tlai,
    tsai,
    frac_sno,
    obu,
    fclay,
    z0a,
    frac_rock,
    frac_veg,
    Kc_map,
    is_soil,
    ovr_src_snk_mss=OVR_SRC_SNK_MSS,
):
    """
    Compute the Leung et al. (2023/2024) dust emission flux.

    Implements the full improved dust emission scheme described in
    Leung et al. (2023) and evaluated in Leung et al. (2024) as a
    vectorized NumPy function over a 2-D [x, y] spatial grid.

    Physics pipeline
    ----------------
    1.  Mobilizable bare-ground fraction (fbare) from VAI and snow cover.
        Uses extended VAIthr = 1.0 to capture marginal semiarid sources.
    2.  Bulk density and gravimetric soil water content.
    3.  Soil-moisture threshold factor fm (Fécan et al. 1999; a = 2).
    4.  Liquid fraction (inhibits emission over frozen soil).
    5.  Dry fluid threshold u*ft0 via Shao & Lu (2000) at Dp = 127 µm.
    6.  Wet fluid threshold u*ft = u*ft0 * fm.
    7.  Standardised threshold u*st = u*ft * sqrt(rho_air/rho_air0).
    8.  Impact (dynamic) threshold u*it = Bit * u*ft0_dry (Comola 2019).
    9.  Drag partition (rock + vegetation) → effective u*s = fv * Feff.
    10. Wind-speed std at saltation height (Panofsky 1977).
    11. Intermittency factor eta in [0, 1] (Comola et al. 2019; CTSM form).
    12. Soil erodibility coefficient Cd(u*st) (Kok et al. 2014b).
    13. Fragmentation exponent kappa(u*st), capped at KAPPA_MAX.
    14. Dust emission flux Fd using K14 equation with u*it threshold.
    15. Apply upscaling correction map Kc (Leung et al. 2024 Eq. 12).
    16. Partition total flux into four transport bins via lognormal overlaps.

    Parameters
    ----------
    fv           : np.ndarray [x, y]   friction velocity [m/s]
    forc_rho     : np.ndarray [x, y]   near-surface air density [kg/m³]
    h2osoi_vol   : np.ndarray [x, y]   volumetric soil moisture [m³/m³]
    h2osoi_liq   : np.ndarray [x, y]   liquid soil water [kg/m²]
    h2osoi_ice   : np.ndarray [x, y]   frozen soil water [kg/m²]
    watsat       : np.ndarray [x, y]   saturated volumetric soil water [m³/m³]
    tlai         : np.ndarray [x, y]   one-sided leaf area index [-]
    tsai         : np.ndarray [x, y]   one-sided stem area index [-]
    frac_sno     : np.ndarray [x, y]   snow cover fraction [-], in [0, 1]
    obu          : np.ndarray [x, y]   Obukhov length [m]
    fclay        : np.ndarray [x, y]   clay mass fraction [-] (raw; used for
                                       gwc_thr, and via clay_mass_fraction_leung
                                       for the emission equation)
    z0a          : np.ndarray [x, y]   aeolian roughness length [m]
    frac_rock    : np.ndarray [x, y]   fractional rock area Ar [-]
    frac_veg     : np.ndarray [x, y]   fractional vegetation area Av [-]
    Kc_map       : np.ndarray [x, y]   upscaling correction factor (Leung 2024)
    is_soil      : np.ndarray [x, y]   bool, True for soil/crop land type
    ovr_src_snk_mss : np.ndarray [3,4] source-to-sink mass overlap fractions

    Returns
    -------
    F_p : np.ndarray [x, y, ndst]   dust emission flux [kg m-2 s-1]
                                     ndst = 4 bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm
    """
    nx, ny = fv.shape
    # Clay term in the emission equation: CTSM MassFracClayLeung2023
    #   fclay' = 0.1 + 0.5 * min(fclay, 0.2)   (in [0.1, 0.2])
    fclay_eff = clay_mass_fraction_leung(fclay)

    # ------------------------------------------------------------------
    # Step 1: Mobilizable bare-ground fraction (VAIthr = 1.0)
    # ------------------------------------------------------------------
    vai = tlai + tsai
    fbare = np.where(
        is_soil,
        np.maximum(0.0, 1.0 - vai / VAI_MBL_THR) * (1.0 - frac_sno),
        0.0,
    )

    # ------------------------------------------------------------------
    # Step 2-3: Soil moisture factor fm
    # ------------------------------------------------------------------
    fm, gwc_sfc, gwc_thr = moisture_threshold_factor(h2osoi_vol, watsat, fclay)

    # ------------------------------------------------------------------
    # Step 4: Liquid fraction (inhibits emission over frozen soil)
    # ------------------------------------------------------------------
    liqfrac = np.clip(
        h2osoi_liq / (h2osoi_ice + h2osoi_liq + 1.0e-6),
        0.0,
        1.0,
    )

    # ------------------------------------------------------------------
    # Step 5-7: Thresholds
    # ------------------------------------------------------------------
    u_ft0_dry = sl2000_u_ft0_dry(forc_rho)  # dry fluid threshold [m/s]
    u_ft = u_ft0_dry * fm  # wet fluid threshold [m/s]
    u_st = u_ft * np.sqrt(forc_rho / RHO_AIR0)  # standardised threshold [m/s]
    u_it = BIT * u_ft0_dry  # impact threshold [m/s]

    # ------------------------------------------------------------------
    # Step 8-9: Drag partition → effective soil-surface u*s
    # ------------------------------------------------------------------
    Feff, feff_r, feff_v = drag_partition(z0a, vai, frac_rock, frac_veg)
    u_s = fv * Feff  # effective u*s [m/s]

    # ------------------------------------------------------------------
    # Step 10-11: Wind-speed fluctuations and intermittency factor
    # ------------------------------------------------------------------
    sigma_U = turbulent_sigma_us(u_s, obu)  # wind-speed std at Z_SAL [m/s]
    eta = intermittency_factor(u_s, sigma_U, u_it, u_ft)

    # ------------------------------------------------------------------
    # Step 12: Soil erodibility coefficient Cd (Kok et al. 2014b Eq. 5b)
    # ------------------------------------------------------------------
    u_st_safe = np.maximum(u_st, 1.0e-6)
    Cd = CD0 * np.exp(-CE * (u_st_safe - U_ST0) / U_ST0)

    # ------------------------------------------------------------------
    # Step 13: Fragmentation exponent kappa (capped at KAPPA_MAX)
    # ------------------------------------------------------------------
    kappa_raw = CK * (u_st_safe - U_ST0) / U_ST0
    kappa = np.minimum(kappa_raw, KAPPA_MAX)

    # ------------------------------------------------------------------
    # Step 14: Dust emission flux (Leung et al. 2023; CTSM DustEmisLeung2023)
    #   Fd = eta * Ctune * Cd * fbare * fclay' * liqfrac
    #        * rho_air * (u*s^2 - u*it^2) / u*it * (u*s/u*it)^kappa
    #   only where u*s > u*it;  fclay' = 0.1 + 0.5*min(fclay, 0.2)  (CTSM)
    #   Units: [kg m-3][m2 s-2]/[m s-1] = [kg m-2 s-1]  (no 1/g factor)
    # ------------------------------------------------------------------
    emit_mask = (fbare > 0.0) & (u_s > u_it)
    u_s_safe = np.maximum(u_s, 1.0e-6)
    u_it_safe = np.maximum(u_it, 1.0e-6)

    flux_term = forc_rho * (u_s_safe**2 - u_it_safe**2) / u_it_safe  # [kg m-2 s-1]
    shape_term = (u_s_safe / u_it_safe) ** kappa  # [-]

    Fd = np.where(
        emit_mask,
        eta * CTUNE * Cd * fbare * fclay_eff * liqfrac * flux_term * shape_term,
        0.0,
    )
    Fd = np.maximum(Fd, 0.0)  # guard against numerics

    # ------------------------------------------------------------------
    # Step 15: Upscaling correction map (Leung et al. 2024, Eq. 12)
    # ------------------------------------------------------------------
    Fd_corrected = Fd * np.maximum(Kc_map, 0.0)

    # ------------------------------------------------------------------
    # Step 16: Partition total flux into four transport bins
    # ------------------------------------------------------------------
    F_p = np.zeros((nx, ny, NDST))
    for n in range(NDST):
        bin_frac = ovr_src_snk_mss[:, n].sum()  # sum over source modes
        F_p[:, :, n] = Fd_corrected * bin_frac

    return F_p  # [kg m-2 s-1], shape [x, y, ndst=4]
