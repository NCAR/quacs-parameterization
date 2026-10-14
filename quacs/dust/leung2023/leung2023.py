"""
Dust Emission Scheme: Leung et al. (2023) / Kok et al. (2014) process-based, scale-aware
=========================================================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Leung, D. M., et al. (2023). A new process-based and scale-aware desert
    dust emission scheme for global climate models - Part I: Description and
    evaluation against inverse modeling emissions. Atmos. Chem. Phys., 23,
    6487-6523. https://doi.org/10.5194/acp-23-6487-2023   (scheme of record)

    Leung, D. M., et al. (2024). ... Part II: Evaluation in CESM2.
    Atmos. Chem. Phys., 24, 2287-2318.
    https://doi.org/10.5194/acp-24-2287-2024              (evaluation)

    Kok, J. F., et al. (2014). An improved dust emission model - Part 1.
    Atmos. Chem. Phys., 14, 13023-13041.                 (K14 base flux)
    Shao, Y. and Lu, H. (2000). A simple expression for wind erosion
    threshold friction velocity. JGR, 105, 22437-22443.  (S&L00 threshold)
    Comola, F., et al. (2019). Intermittent dust emission ... GRL, 46.
    Okin, G. S. (2008). JGR, 113, F02S10.  Pierre, C., et al. (2014). drag part.
    Fecan, F., et al. (1999). Ann. Geophys., 17, 149-157.  (moisture)

    Fortran reference (CLM5 / CTSM, CESM2.1 fork dustmodelcesm2023):
      DUSTMod.F90  subroutine DustEmission, ACTIVE (uncommented) path only:
        L505-L511  lnd_frc_mbl (bare/mobile fraction)
        L571-L577  Fecan moisture factor
        L593-L596  Shao & Lu (2000) threshold; impact threshold
        L605-L612  standardized threshold, emission coeff Cd, frag exponent
        L623-L641  Okin/Pierre hybrid drag partition on the wind
        L724,L752  K14 impact-threshold vertical flux
        L755       bare-fraction / liqfrac / C_tune scaling
        L778-L811  Comola (2019) intermittency factor
      SoilStateInitTimeConstMod.F90
        L675       gwc_thr = 0.01*2*(0.17*%clay + 0.0014*%clay^2)  (a=2)
        L676       mss_frc_cly_vld = min(%clay/100, 0.20)
      SoilStateInitTimeConstMod.F90 L329/L343-345  static drag fields
        roughfct (F_eff rock), rockfrc (A_r), vegefrc (A_v), sparfrc (A_s)

QUACS Contract
--------------
Standalone module. Dynamic met in, surface flux out, NumPy-vectorized.
Runs as a script (python dust_leung2023.py) or pytest (pytest dust_leung2023.py).
One dependency beyond NumPy: scipy.special.erf, required by the intermittency
factor (Comola 2019). This is the only scheme in the set that needs it.

Physics summary (active K14 impact-threshold path):
    u*ft0  = sqrt(A*(rho_p*g*Dp + gamma/Dp)) / sqrt(rho_a)        [S&L00]
    fm     = sqrt(1 + 1.21*(100*(w - w_t))^0.68)   for w > w_t    [Fecan, a=2]
    u*ft   = u*ft0 * fm                              (wet fluid threshold)
    u*it   = B_it * u*ft0                            (impact threshold; uses DRY)
    u*st   = u*ft * sqrt(rho_a / rho_a0)             (standardized threshold)
    Cd     = Cd0 * exp(-Ce*(u*st - u*st0)/u*st0)     (emission coefficient)
    n      = min(Ca*(u*st - u*st0)/u*st0, 3)         (fragmentation exponent)
    u*s    = u* * F_eff                              (drag partition on the wind)
    F      = Cd * f_clay * rho_a * (u*s^2 - u*it^2)/u*st * (u*s/u*it)^n
                                                     for u*s > u*it
    F      = F * lnd_frc_mbl * C_tune * liqfrac * eta

Three structural points worth knowing before you trust the numbers
(all confirmed against the Fortran, not the paper):
  1. The impact threshold u*it uses the DRY fluid threshold (DUSTMod L596),
     so soil moisture raises u*ft but never the emission gate u*it. Moisture
     enters the flux only through u*ft -> u*st -> (Cd, n).
  2. Dp is 130 um in the code (DUSTMod L593), not the 127 um the papers report.
     Pinned here to the code value.
  3. Atmospheric stability is hardcoded to 0 (DUSTMod L788), so the wind
     standard deviation collapses to u_sd = u*s * 12^0.333 and the
     intermittency factor becomes a closed form in (u*s, u*ft, u*it) alone.
     No Obukhov length or external wind-fluctuation input is required.
  The CTSM driver applies the coastal factor and the pinned four-bin overlap
  table; its total is the sum of the four bins. The pre-split flux is diagnostic.
  Note: cst_slt = 2.61 is declared in DUSTMod (L382) but is dead in this path;
  it appears only in the commented-out Z03 block. Not used here.

Inputs (dynamic, per time step from MPAS-A):
    ustar    : np.ndarray [x, y]   friction velocity u* [m/s]   (driver)
    rhoa     : np.ndarray [x, y]   near-surface air density [kg/m^3]
    gwc      : np.ndarray [x, y]   gravimetric soil moisture [kg/kg]
    lai      : np.ndarray [x, y]   LAI + SAI (vegetation area index) [m2/m2]
    frac_sno : np.ndarray [x, y]   snow-covered ground fraction [-], in [0,1]
    liqfrac  : np.ndarray [x, y]   liquid fraction of soil water [-] (1=unfrozen)

Inputs (static, loaded once at initialization):
    pct_clay : np.ndarray [x, y]   clay percentage [0-100]; derives f_clay & w_t
    rockfrc  : np.ndarray [x, y]   rock area fraction A_r [-]
    vegefrc  : np.ndarray [x, y]   shrub/crop area fraction A_v [-]
    sparfrc  : np.ndarray [x, y]   sparse-vegetation area fraction A_s [-]
    roughfct : np.ndarray [x, y]   rock drag partition factor f_eff,r [-] (F_eff)

Output:
    dust_emission returns pre-split surface flux. leung_dust_driver returns
    four bin fluxes and their sum after CTSM soil and patch input steps.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field, fields

import numpy as np
from scipy.special import erf

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80616  # m/s^2  (CLM clm_varcon grav)
RDAIR = 6.02214e26 * 1.38065e-23 / 28.966

# ---------------------------------------------------------------------------
# Scheme constants -- pinned to DUSTMod.F90 unless noted
# ---------------------------------------------------------------------------
# Shao & Lu (2000) threshold (DUSTMod L593)
SL00_A = 0.0123  # [-]       interparticle-force coefficient A
SL00_GAM = 1.65e-4  # [kg/s^2]  interparticle-force constant gamma
DP = 130.0e-6  # [m]       median soil diameter (CODE value; papers say 127 um)
DNS_SLT = 2650.0  # [kg/m^3]  saltator density (DUSTMod L393)

# Fecan et al. (1999) moisture (SoilStateInitTimeConstMod L675, a = 2)
FECAN_A = 2.0  # [-]       moisture tunable a

# Kok et al. (2014) coefficients (DUSTMod L387-L392)
CD0 = 4.4e-5  # [-]
CE = 2.0  # [-]
CA = 2.7  # [-]
C_TUNE = 0.05  # [-]       global vertical-flux tuning constant
UST0_MIN = 0.16  # [m/s]     minimum standardized threshold
FORC_RHO_STD = 1.2250  # [kg/m^3]  reference air density

# Intermittency / impact threshold (DUSTMod L386, L395-L399)
VAI_MBL_THR = 1.0  # [m2/m2]   VAI above which the cell is non-mobile
B_IT = 0.82  # [-]       u*it / u*ft0
KARMAN = 0.4  # [-]       von Karman constant
U_0 = 0.32  # [-]       Okin SSR in the immediate lee of a plant
C_E = 4.8  # [-]       Okin velocity-recovery e-folding distance
Z_SAL = 0.1  # [m]       saltation height
Z0_SAL = 1.0e-4  # [m]       assumed aeolian roughness length at saltation

# Literal numeric constants used verbatim by the Fortran. Kept as-is (not 1/3,
# not sqrt(2)) so a scalar transcription reproduces the Fortran to machine
# precision -- DUSTMod L790 uses 0.333, L806/L808 use 1.414.
SD_EXP = 0.333  # DUSTMod L790  (12 - 0.5*stblty)**0.333, stblty := 0
SQRT2_LIT = 1.414  # DUSTMod L806/L808  sqrt(2) literal
CBRT_LIT = 0.3333  # DUSTMod L633  drag-partition cube root written as **0.3333
# (a Fortran approximation to 1/3; kept verbatim so the
#  transcription matches CLM bit-for-bit near threshold)


def _setting(value, source, what):
    return field(default=value, metadata={"source": source, "what": what})


@dataclass(frozen=True)
class LeungSettings:
    sl00_a: float = _setting(SL00_A, "DUSTMod.F90:593", "dry threshold scale")
    sl00_gam: float = _setting(SL00_GAM, "DUSTMod.F90:593", "interparticle force")
    dp: float = _setting(DP, "DUSTMod.F90:593", "soil particle diameter")
    dns_slt: float = _setting(DNS_SLT, "DUSTMod.F90:393", "saltator density")
    fecan_a: float = _setting(
        FECAN_A, "SoilStateInitTimeConstMod.F90:675", "moisture threshold scale"
    )
    cd0: float = _setting(CD0, "DUSTMod.F90:387", "emission scale")
    ce: float = _setting(CE, "DUSTMod.F90:389", "emission threshold slope")
    ca: float = _setting(CA, "DUSTMod.F90:388", "fragmentation slope")
    c_tune: float = _setting(C_TUNE, "DUSTMod.F90:390", "vertical flux scale")
    ust0_min: float = _setting(UST0_MIN, "DUSTMod.F90:391", "standard threshold reference")
    forc_rho_std: float = _setting(FORC_RHO_STD, "DUSTMod.F90:392", "standard air density")
    vai_mbl_thr: float = _setting(VAI_MBL_THR, "DUSTMod.F90:386", "vegetation mobility cutoff")
    b_it: float = _setting(B_IT, "DUSTMod.F90:395", "impact threshold ratio")
    u_0: float = _setting(U_0, "DUSTMod.F90:398", "plant lee stress ratio")
    c_e: float = _setting(C_E, "DUSTMod.F90:399", "plant lee recovery")
    z_sal: float = _setting(Z_SAL, "DUSTMod.F90:778", "saltation height")
    z0_sal: float = _setting(Z0_SAL, "DUSTMod.F90:778", "saltation roughness")
    sd_exp: float = _setting(SD_EXP, "DUSTMod.F90:790", "wind variance exponent")
    sqrt2_lit: float = _setting(SQRT2_LIT, "DUSTMod.F90:806", "erf argument scale")
    cbrt_lit: float = _setting(CBRT_LIT, "DUSTMod.F90:633", "drag cube root exponent")
    clay_linear: float = _setting(0.17, "SoilStateInitTimeConstMod.F90:675", "clay threshold slope")
    clay_quadratic: float = _setting(
        0.0014, "SoilStateInitTimeConstMod.F90:675", "clay threshold curvature"
    )
    wet_scale: float = _setting(1.21, "DUSTMod.F90:574", "wet threshold scale")
    wet_exp: float = _setting(0.68, "DUSTMod.F90:574", "wet threshold exponent")
    lai_offset: float = _setting(0.01, "DUSTMod.F90:625", "drag LAI offset")
    lai_cap: float = _setting(1.0, "DUSTMod.F90:627", "drag LAI ceiling")
    gap_multiplier: float = _setting(2.0, "DUSTMod.F90:630", "gap length scale")
    drag_gate: float = _setting(1.0, "DUSTMod.F90:623", "drag VAI cutoff")
    variance_neutral: float = _setting(12.0, "DUSTMod.F90:789", "neutral wind variance")
    stability_coeff: float = _setting(0.5, "DUSTMod.F90:789", "stability wind variance slope")
    sd_floor: float = _setting(0.001, "DUSTMod.F90:791", "wind variance branch threshold")
    clay_cap: float = _setting(0.20, "SoilStateInitTimeConstMod.F90:676", "valid clay cap")
    frag_cap: float = _setting(3.0, "DUSTMod.F90:610", "fragmentation exponent cap")
    coast_cutoff: float = _setting(0.5, "DUSTMod.F90:758", "coastal land cutoff")
    coast_scale: float = _setting(1e-3, "DUSTMod.F90:759", "coastal flux scale")
    soil_particle_density: float = _setting(2700.0, "DUSTMod.F90:571", "soil particle density")
    liquid_guard: float = _setting(1e-6, "DUSTMod.F90:583", "liquid fraction guard")
    weight_tolerance: float = _setting(1e-6, "DUSTMod.F90:475", "patch weight tolerance")
    organic_max: float = _setting(
        130.0,
        "SoilStateInitTimeConstMod.F90:240; clm5_params.c171117.nc",
        "maximum organic density",
    )
    zsapric: float = _setting(0.5, "SoilStateInitTimeConstMod.F90:133", "sapric depth")
    zsoi1: float = _setting(0.01, "initVerticalMod.F90:260-281", "top soil midpoint")
    cosby_base: float = _setting(0.489, "FuncPedotransferMod.F90:98", "mineral porosity intercept")
    cosby_sand: float = _setting(
        0.00126, "FuncPedotransferMod.F90:98", "mineral porosity sand slope"
    )
    organic_porosity_base: float = _setting(
        0.93, "SoilStateInitTimeConstMod.F90:514", "organic porosity intercept"
    )
    organic_porosity_slope: float = _setting(
        0.1, "SoilStateInitTimeConstMod.F90:514", "organic porosity depth slope"
    )
    organic_porosity_min: float = _setting(
        0.83, "SoilStateInitTimeConstMod.F90:514", "organic porosity floor"
    )
    source_medians: tuple = _setting(
        (0.832e-6, 4.82e-6, 19.38e-6), "DUSTMod.F90:1080", "source median diameters"
    )
    source_gsds: tuple = _setting(
        (2.10, 1.90, 1.60), "DUSTMod.F90:1081", "source geometric deviations"
    )
    source_fractions: tuple = _setting(
        (0.036, 0.957, 0.007), "DUSTMod.F90:1082", "source mass fractions"
    )
    bin_grid: tuple = _setting(
        (0.1e-6, 1e-6, 2.5e-6, 5e-6, 10e-6), "DUSTMod.F90:1060-1061", "dust bin boundaries"
    )
    ndst: int = _setting(4, "DUSTMod.F90:1145", "required bin count")
    sz_nbr: int = _setting(200, "clm_varpar.F90:46", "size integration intervals")
    gsd_anl: float = _setting(2.0, "DUSTMod.F90:1152", "deposition size distribution width")
    dns_aer: float = _setting(2500.0, "DUSTMod.F90:1154", "aerosol density")
    dmt_vma: float = _setting(3.5e-6, "DUSTMod.F90:1158", "deposition mass median diameter")
    eps_max: float = _setting(1e-4, "DUSTMod.F90:1214", "settling convergence tolerance")
    setup_t: float = _setting(295.0, "DUSTMod.F90:1215", "setup air temperature")
    setup_p: float = _setting(1e5, "DUSTMod.F90:1215", "setup air pressure")
    mu_scale: float = _setting(1.72e-5, "DUSTMod.F90:925", "air viscosity scale")
    mu_tref: float = _setting(273.0, "DUSTMod.F90:925", "air viscosity temperature reference")
    mu_s: float = _setting(393.0, "DUSTMod.F90:925", "air viscosity multiplier")
    mu_offset: float = _setting(120.0, "DUSTMod.F90:926", "air viscosity temperature offset")
    slip_a: float = _setting(1.257, "DUSTMod.F90:932", "slip correction intercept")
    slip_b: float = _setting(0.4, "DUSTMod.F90:932", "slip correction exponential scale")
    slip_exp: float = _setting(1.1, "DUSTMod.F90:932", "slip correction exponent")
    schmidt_exp: float = _setting(-2 / 3, "DUSTMod.F90:895", "land Schmidt exponent")
    impaction_exp: float = _setting(3.0, "DUSTMod.F90:961", "impaction Stokes exponent")
    re_cut1: float = _setting(0.1, "DUSTMod.F90:1257", "first Reynolds cutoff")
    re_cut2: float = _setting(2.0, "DUSTMod.F90:1259", "second Reynolds cutoff")
    re_cut3: float = _setting(500.0, "DUSTMod.F90:1264", "third Reynolds cutoff")
    re_max: float = _setting(2e5, "DUSTMod.F90:1266", "maximum Reynolds number")
    drag_low_linear: float = _setting(3 / 16, "DUSTMod.F90:1261", "low Reynolds drag slope")
    drag_low_log: float = _setting(9 / 160, "DUSTMod.F90:1263", "low Reynolds log term")
    drag_mid_scale: float = _setting(0.15, "DUSTMod.F90:1265", "mid Reynolds drag scale")
    drag_mid_exp: float = _setting(0.687, "DUSTMod.F90:1265", "mid Reynolds exponent")
    drag_high: float = _setting(0.44, "DUSTMod.F90:1267", "high Reynolds drag coefficient")
    convergence_average: int = _setting(12, "DUSTMod.F90:1281", "settling averaging iteration")
    convergence_warn: int = _setting(20, "DUSTMod.F90:1286", "settling warning iteration")
    use_drag_partition: bool = _setting(True, "dust_leung2023.py:438", "enable vegetation drag")
    use_intermittency: bool = _setting(True, "dust_leung2023.py:457", "enable wind intermittency")


def tunable_names(settings_class):
    return tuple(item.name for item in fields(settings_class))


class LeungInputWarning(UserWarning):
    """Input outside the physical domain of the pinned CTSM calculation."""


# ---------------------------------------------------------------------------
# Threshold: Shao & Lu (2000)  (DUSTMod L593-L596)
# ---------------------------------------------------------------------------
def dry_fluid_threshold(rhoa, settings=None):
    """
    Dry fluid threshold friction velocity u*ft0 [m/s], Shao & Lu (2000).

        u*ft0 = sqrt(A*(rho_p*g*Dp + gamma/Dp)) / sqrt(rho_a)

    Parameters
    ----------
    rhoa : np.ndarray [x, y]   air density [kg/m^3]

    Returns
    -------
    u_ft0 : np.ndarray [x, y]  dry fluid threshold [m/s]
    """
    s = LeungSettings() if settings is None else settings
    num = np.sqrt(s.sl00_a * (s.dns_slt * GRAV * s.dp + s.sl00_gam / s.dp))
    return num / np.sqrt(rhoa)


# ---------------------------------------------------------------------------
# Moisture factor: Fecan et al. (1999)  (DUSTMod L571-L577; threshold L675)
# ---------------------------------------------------------------------------
def moisture_factor(gwc, pct_clay, settings=None):
    """
    Fecan moisture enhancement factor fm (>= 1) and the threshold w_t.

        w_t = 0.01 * a * (0.17*%clay + 0.0014*%clay^2)
        fm  = sqrt(1 + 1.21*(100*(w - w_t))^0.68)   for w > w_t, else 1

    Note the sqrt form matches the original Fecan (1999) and the CLM code;
    paper Eq. (2a) is written without the sqrt. Code is the source of record.
    """
    s = LeungSettings() if settings is None else settings
    w_t = 0.01 * s.fecan_a * (s.clay_linear * pct_clay + s.clay_quadratic * pct_clay * pct_clay)
    excess = np.maximum(100.0 * (gwc - w_t), 0.0)
    fm = np.where(gwc > w_t, np.sqrt(1.0 + s.wet_scale * excess**s.wet_exp), 1.0)
    return fm, w_t


# ---------------------------------------------------------------------------
# Hybrid drag partition: Okin (2008) + Pierre (2014) + rocks  (DUSTMod L623-L641)
# ---------------------------------------------------------------------------
def drag_partition(lai, rockfrc, vegefrc, sparfrc, roughfct, settings=None):
    """
    Hybrid drag partition factor F_eff [-] combining rock roughness (precomputed
    roughfct = f_eff,r) and vegetation (Okin SSR via Pierre's gap-length form).

        K     = 2*(1/lai_c - 1),  lai_c = min(lai + 0.01, 1)
        ssr   = (K + U_0*c_e) / (K + c_e)
        F_eff = (A_r*f_eff,r^3 + (A_v + A_s)*ssr^3)^(1/3)

    Applied only where the cell is mobile and lai <= 1 (caller enforces the
    mobility gate); elsewhere the wind is used undrag-partitioned (F_eff = 1).
    """
    s = LeungSettings() if settings is None else settings
    lai_c = np.minimum(lai + s.lai_offset, s.lai_cap)
    K = s.gap_multiplier * (1.0 / lai_c - 1.0)
    ssr = (K + s.u_0 * s.c_e) / (K + s.c_e)
    f_eff = (rockfrc * roughfct**3 + (vegefrc + sparfrc) * ssr**3) ** s.cbrt_lit
    return f_eff


# ---------------------------------------------------------------------------
# Intermittency factor: Comola et al. (2019)  (DUSTMod L778-L811, stblty := 0)
# ---------------------------------------------------------------------------
def intermittency_factor(u_s, u_ft, u_it, settings=None, stability=0.0):
    """
    Intermittency factor eta in [0, 1], the time fraction of active emission.

    With stability hardcoded to 0 (DUSTMod L788) this is a closed form in the
    soil-surface friction velocity and the two thresholds, all translated to
    the saltation height z = 0.1 m by a neutral log law.

        u_bar = (u_s /k) * ln(z_sal/z0_sal)
        u_sd  =  u_s     * 12^0.333
        u_ftz = (u_ft/k) * ln(z_sal/z0_sal)
        u_itz = (u_it/k) * ln(z_sal/z0_sal)
        R     = 1 / (exp((u_ftz^2 - u_itz^2 - 2*u_bar*(u_ftz - u_itz))
                          / (2*u_sd^2)) + 1)
        P_ft  = 0.5*(1 + erf((u_ftz - u_bar)/(1.414*u_sd)))
        P_it  = 0.5*(1 + erf((u_itz - u_bar)/(1.414*u_sd)))
        eta   = 1 - P_ft + R*(P_ft - P_it)
    """
    s = LeungSettings() if settings is None else settings
    log_ratio = np.log(s.z_sal / s.z0_sal)
    u_bar = (u_s / KARMAN) * log_ratio
    variance = s.variance_neutral - s.stability_coeff * stability
    u_sd = u_s * variance**s.sd_exp if variance >= s.sd_floor else np.full_like(u_s, s.sd_floor)
    u_ftz = (u_ft / KARMAN) * log_ratio
    u_itz = (u_it / KARMAN) * log_ratio

    # exp may overflow to +inf for very small u_sd; the limit R -> 0 is correct,
    # so suppress the benign overflow rather than let it warn in SA sweeps.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        R = 1.0 / (
            np.exp((u_ftz**2 - u_itz**2 - 2.0 * u_bar * (u_ftz - u_itz)) / (2.0 * u_sd**2)) + 1.0
        )
        P_ft = 0.5 * (1.0 + erf((u_ftz - u_bar) / (s.sqrt2_lit * u_sd)))
        P_it = 0.5 * (1.0 + erf((u_itz - u_bar) / (s.sqrt2_lit * u_sd)))
    eta = 1.0 - P_ft + R * (P_ft - P_it)
    return eta


# ---------------------------------------------------------------------------
# Initialization: load static drag-partition and soil fields once at start
# ---------------------------------------------------------------------------
def initialize(clay_file, rockfrc_file, vegefrc_file, sparfrc_file, roughfct_file):
    """
    Load and validate static input fields. Called once before the time loop.

    Returns
    -------
    pct_clay, rockfrc, vegefrc, sparfrc, roughfct : np.ndarray [x, y]
    """
    files = [clay_file, rockfrc_file, vegefrc_file, sparfrc_file, roughfct_file]
    labels = ["clay_file", "rockfrc_file", "vegefrc_file", "sparfrc_file", "roughfct_file"]
    for path, label in zip(files, labels):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_leung2023] {label} not found: {path}")

    pct_clay = np.load(clay_file)
    rockfrc = np.load(rockfrc_file)
    vegefrc = np.load(vegefrc_file)
    sparfrc = np.load(sparfrc_file)
    roughfct = np.load(roughfct_file)

    shp = pct_clay.shape
    if pct_clay.ndim != 2:
        raise ValueError(f"[dust_leung2023] pct_clay must be 2D [x,y], got {shp}")
    if not (pct_clay.min() >= 0.0 and pct_clay.max() <= 100.0):
        raise ValueError(
            f"[dust_leung2023] pct_clay must be a percentage in [0,100]; "
            f"got [{pct_clay.min():.2f}, {pct_clay.max():.2f}]"
        )
    for arr, label in [
        (rockfrc, "rockfrc"),
        (vegefrc, "vegefrc"),
        (sparfrc, "sparfrc"),
        (roughfct, "roughfct"),
    ]:
        if arr.shape != shp:
            raise ValueError(f"[dust_leung2023] {label} shape {arr.shape} != pct_clay shape {shp}")
        if not (arr.min() >= 0.0 and arr.max() <= 1.0):
            raise ValueError(
                f"[dust_leung2023] {label} must be in [0,1]; got [{arr.min():.3f}, {arr.max():.3f}]"
            )

    return pct_clay, rockfrc, vegefrc, sparfrc, roughfct


# ---------------------------------------------------------------------------
# Core parameterization returns the pre-split flux used by the CTSM driver.
# ---------------------------------------------------------------------------
def dust_emission(
    ustar,
    rhoa,
    gwc,
    lai,
    frac_sno,
    pct_clay,
    rockfrc,
    vegefrc,
    sparfrc,
    roughfct,
    liqfrac=1.0,
    settings=None,
    return_diagnostics=False,
    mobile_fraction=None,
):
    """
    Total vertical dust emission flux, Leung (2023) active K14 impact-threshold path.

    Parameters
    ----------
    ustar, rhoa, gwc, lai, frac_sno : np.ndarray [x, y]   dynamic met (see module doc)
    pct_clay, rockfrc, vegefrc, sparfrc, roughfct : np.ndarray [x, y]   static
    liqfrac : np.ndarray [x, y] or float   liquid soil-water fraction [-]
    settings : LeungSettings   source-pinned coefficients and two model switches
    return_diagnostics : bool   also return a dict of intermediate fields

    Returns
    -------
    F : np.ndarray [x, y]   total dust emission flux [kg m^-2 s^-1]
        (plus diagnostics dict if return_diagnostics)
    """
    s = LeungSettings() if settings is None else settings
    ustar = np.asarray(ustar, dtype=float)
    rhoa = np.asarray(rhoa, dtype=float)

    # --- valid clay fraction, capped at 0.20 (SoilState L676) ----------------
    f_clay = np.minimum(pct_clay / 100.0, s.clay_cap)

    # --- bare / mobile fraction (DUSTMod L505-L511) --------------------------
    lnd_frc_mbl = (
        np.where(lai < s.vai_mbl_thr, 1.0 - lai / s.vai_mbl_thr, 0.0) * (1.0 - frac_sno)
        if mobile_fraction is None
        else mobile_fraction
    )
    mobile = lnd_frc_mbl > 0.0

    # --- thresholds: S&L00 dry, Fecan-wet fluid, impact (DUSTMod L593-L596) --
    u_ft0 = dry_fluid_threshold(rhoa, settings=s)
    fm, _ = moisture_factor(gwc, pct_clay, settings=s)
    u_ft = u_ft0 * fm  # wet fluid threshold
    tmp2 = np.sqrt(s.sl00_a * (s.dns_slt * GRAV * s.dp + s.sl00_gam / s.dp))
    u_it = s.b_it * tmp2 / np.sqrt(rhoa)

    # --- standardized threshold, Cd, fragmentation exponent (L605-L612) ------
    u_st = u_ft * np.sqrt(rhoa / s.forc_rho_std)
    Cd = s.cd0 * np.exp(-s.ce * (u_st - s.ust0_min) / s.ust0_min)
    frag_expt = s.ca * (u_st - s.ust0_min) / s.ust0_min
    frag_expt = np.minimum(frag_expt, s.frag_cap)

    # --- drag partition on the wind: u*s = u* * F_eff (L623-L641) ------------
    if s.use_drag_partition:
        f_eff = drag_partition(lai, rockfrc, vegefrc, sparfrc, roughfct, settings=s)
        apply = mobile & (lai <= s.drag_gate)
        u_s = np.where(apply, ustar * f_eff, ustar)
    else:
        u_s = ustar

    # --- K14 impact-threshold flux, gated on u*s > u*it (L724, L752) ---------
    gate = u_s > u_it
    with np.errstate(invalid="ignore", divide="ignore"):
        flux = (
            Cd
            * f_clay
            * rhoa
            * (u_s**2 - u_it**2)
            / u_st
            * np.where(gate, (u_s / u_it), 0.0) ** frag_expt
        )
    flux = np.where(gate, flux, 0.0)

    # --- bare-fraction / tuning / liquid-fraction scaling (L755) -------------
    flux = flux * lnd_frc_mbl * s.c_tune * liqfrac

    # --- intermittency factor (L778-L811); NaN-safe like L814-L818 -----------
    if s.use_intermittency:
        eta = intermittency_factor(u_s, u_ft, u_it, settings=s)
        flux = np.where(np.isnan(eta), flux, flux * eta)
    else:
        eta = np.ones_like(flux)

    F = flux

    if return_diagnostics:
        return F, dict(
            u_ft0=u_ft0,
            u_ft=u_ft,
            u_it=u_it,
            u_st=u_st,
            Cd=Cd,
            frag_expt=frag_expt,
            u_s=u_s,
            lnd_frc_mbl=lnd_frc_mbl,
            fm=fm,
            eta=eta,
            f_clay=f_clay,
        )
    return F


def init_dust_vars(settings=None):
    """Pinned InitDustVars size overlap, mean diameters, and Stokes correction."""
    s = LeungSettings() if settings is None else settings
    if s.ndst != 4 or len(s.bin_grid) != s.ndst + 1:
        raise ValueError("DUSTMod.F90:1145 ndst must equal 4")
    overlap = np.empty((len(s.source_medians), s.ndst))
    for m, (median, gsd, fraction) in enumerate(
        zip(s.source_medians, s.source_gsds, s.source_fractions)
    ):
        scale = np.sqrt(2.0) * np.log(gsd)
        for n in range(s.ndst):
            overlap[m, n] = (
                0.5
                * (
                    erf(np.log(s.bin_grid[n + 1] / median) / scale)
                    - erf(np.log(s.bin_grid[n] / median) / scale)
                )
                * fraction
            )

    ln_gsd = np.log(s.gsd_anl)
    dmt_nma = s.dmt_vma * np.exp(-3.0 * ln_gsd * ln_gsd)
    dmt_vwr = np.empty(s.ndst)
    for n in range(s.ndst):
        ratio = (s.bin_grid[n + 1] / s.bin_grid[n]) ** (1.0 / s.sz_nbr)
        sz_min = s.bin_grid[n]
        weighted = 0.0
        volume = 0.0
        for m in range(s.sz_nbr):
            sz_max = sz_min * ratio if m < s.sz_nbr - 1 else s.bin_grid[n + 1]
            center = 0.5 * (sz_min + sz_max)
            tmp = np.log(center / dmt_nma) / ln_gsd
            density = np.exp(-0.5 * tmp * tmp) / center / (ln_gsd * np.sqrt(2.0 * np.pi))
            part = np.pi / 6.0 * center**3.0 * density * (sz_max - sz_min)
            weighted += center * part
            volume += part
            sz_min = sz_max
        dmt_vwr[n] = weighted / volume

    rho = s.setup_p / (s.setup_t * RDAIR)
    mu = s.mu_scale * (s.setup_t / s.mu_tref) ** 1.5 * s.mu_s / (s.setup_t + s.mu_offset)
    mfp = 2.0 * mu / (s.setup_p * np.sqrt(8.0 / (np.pi * RDAIR * s.setup_t)))
    nu = mu / rho
    slip = (
        1.0
        + 2.0 * mfp * (s.slip_a + s.slip_b * np.exp(-s.slip_exp * dmt_vwr / (2.0 * mfp))) / dmt_vwr
    )
    stokes = (1.0 / 18.0) * dmt_vwr**2 * s.dns_aer * GRAV * slip / mu
    correction = np.empty(s.ndst)
    for m in range(s.ndst):
        terminal = stokes[m]
        for iteration in range(max(s.convergence_warn + 3, 25)):
            old = terminal
            reynolds = old * dmt_vwr[m] / nu
            if reynolds < s.re_cut1:
                drag = 24.0 / reynolds
            elif reynolds < s.re_cut2:
                drag = (
                    24.0
                    / reynolds
                    * (
                        1.0
                        + s.drag_low_linear * reynolds
                        + s.drag_low_log * reynolds * reynolds * np.log(2.0 * reynolds)
                    )
                )
            elif reynolds < s.re_cut3:
                drag = 24.0 / reynolds * (1.0 + s.drag_mid_scale * reynolds**s.drag_mid_exp)
            elif reynolds < s.re_max:
                drag = s.drag_high
            else:
                raise ValueError("DUSTMod.F90:1270 Reynolds number too large")
            terminal = np.sqrt(4.0 * GRAV * dmt_vwr[m] * slip[m] * s.dns_aer / (3.0 * drag * rho))
            error = abs((terminal - old) / terminal)
            if iteration == s.convergence_average:
                terminal = 0.5 * (terminal + old)
            if iteration > s.convergence_warn:
                warnings.warn(
                    "DUSTMod.F90:1286 terminal velocity not converging", LeungInputWarning
                )
                break
            if error <= s.eps_max:
                break
        correction[m] = terminal / stokes[m]
    return {
        "ovr_src_snk_mss": overlap,
        "dmt_vwr": dmt_vwr,
        "stk_crc": correction,
        "dns_aer": s.dns_aer,
    }


def dust_dry_dep(ustar, rhoa, forc_pbot, forc_t, ram1, settings=None, setup=None):
    """Pinned DustDryDep turbulent deposition velocity for four dust bins."""
    s = LeungSettings() if settings is None else settings
    setup = init_dust_vars(s) if setup is None else setup
    diameter = setup["dmt_vwr"]
    mu = s.mu_scale * (forc_t / s.mu_tref) ** 1.5 * s.mu_s / (forc_t + s.mu_offset)
    mfp = 2.0 * mu / (forc_pbot * np.sqrt(8.0 / (np.pi * RDAIR * forc_t)))
    nu = mu / rhoa
    slip = (
        1.0
        + 2.0
        * mfp[..., None]
        * (s.slip_a + s.slip_b * np.exp(-s.slip_exp * diameter / (2.0 * mfp[..., None])))
        / diameter
    )
    grav_velocity = (
        (1.0 / 18.0) * diameter**2 * s.dns_aer * GRAV * slip / mu[..., None] * setup["stk_crc"]
    )
    stokes_number = grav_velocity * ustar[..., None] * ustar[..., None] / (GRAV * nu[..., None])
    diffusivity = 1.38065e-23 * forc_t[..., None] * slip / (3.0 * np.pi * mu[..., None] * diameter)
    schmidt = nu[..., None] / diffusivity
    with np.errstate(divide="ignore", invalid="ignore"):
        term = schmidt**s.schmidt_exp + 10.0 ** (-s.impaction_exp / stokes_number)
        resistance = 1.0 / (term * ustar[..., None])
        turbulent_resistance = (
            ram1[..., None] + resistance + ram1[..., None] * resistance * grav_velocity
        )
        velocity = 1.0 / turbulent_resistance
    return {"vlc_trb": velocity, "vlc_grv": grav_velocity, "slp_crc": slip, "rss_lmn": resistance}


def leung_dust_driver(
    ustar,
    rhoa,
    h2osoi_vol,
    tlai,
    tsai,
    frac_sno,
    pct_sand,
    pct_clay,
    organic,
    h2osoi_liq,
    h2osoi_ice,
    land_fraction,
    landunit_type,
    patch_active,
    patch_wtlunit,
    rockfrc,
    vegefrc,
    sparfrc,
    roughfct,
    forc_pbot,
    forc_t,
    ram1,
    settings=None,
    lake=False,
    landunit_id=None,
    pct_sand_layers=None,
    pct_clay_layers=None,
    setup=None,
    initial_value=0.0,
):
    """CTSM soil and patch steps, followed by the shared flux core and four bins."""
    s = LeungSettings() if settings is None else settings
    shape = np.shape(ustar)
    sand = np.broadcast_to(np.asarray(pct_sand, dtype=float), shape).copy()
    clay = np.broadcast_to(np.asarray(pct_clay, dtype=float), shape).copy()
    zero_surface = sand + clay == 0.0
    if pct_sand_layers is not None and pct_clay_layers is not None:
        layers = np.asarray(pct_sand_layers) + np.asarray(pct_clay_layers)
        if np.any(zero_surface & np.any(layers != 0.0, axis=-1)):
            raise ValueError("SoilStateInitTimeConstMod.F90:278 depth points do not sum to zero")
        if np.any(np.any(layers == 0.0, axis=-1) & ~zero_surface):
            raise ValueError("SoilStateInitTimeConstMod.F90:284 depth points sum to zero")
    sand[zero_surface] = 1.0
    clay[zero_surface] = 1.0
    if np.any(sand + clay == 0.0):
        raise ValueError("SoilStateInitTimeConstMod.F90:284 sand plus clay is zero")
    om_frac = np.asarray(organic) / s.organic_max
    om_porosity = max(
        s.organic_porosity_base - s.organic_porosity_slope * s.zsoi1 / s.zsapric,
        s.organic_porosity_min,
    )
    watsat = (1.0 - om_frac) * (s.cosby_base - s.cosby_sand * sand) + om_porosity * om_frac
    gwc = np.asarray(h2osoi_vol) * 1000.0 / ((1.0 - watsat) * s.soil_particle_density)
    liquid_fraction = np.maximum(
        0.0,
        np.minimum(
            1.0,
            np.asarray(h2osoi_liq)
            / (np.asarray(h2osoi_ice) + np.asarray(h2osoi_liq) + s.liquid_guard),
        ),
    )
    active = np.broadcast_to(np.asarray(patch_active, dtype=bool), shape)
    lake = np.broadcast_to(np.asarray(lake, dtype=bool), shape)
    weight = np.broadcast_to(np.asarray(patch_wtlunit), shape)
    ids = (
        np.arange(np.prod(shape)).reshape(shape)
        if landunit_id is None
        else np.broadcast_to(landunit_id, shape)
    )
    vegetation = np.where(lake, 0.0, np.asarray(tlai) + np.asarray(tsai))
    mean_vai = np.empty(shape)
    for unit in np.unique(ids):
        members = ids == unit
        # DUSTMod L455-L472: lake patches keep ttlai = 0 and still count.
        included = members & active & (weight != 0.0)
        total_weight = weight[included].sum()
        if total_weight > 1.0 + s.weight_tolerance:
            raise ValueError("DUSTMod.F90:473 landunit patch weight exceeds one")
        mean_vai[members] = (
            np.sum(vegetation[included] * weight[included]) / total_weight
            if total_weight != 0
            else np.nan
        )
    allowed = (np.asarray(landunit_type) == 1) | (np.asarray(landunit_type) == 2)
    mobile = np.where(
        allowed & ~lake,
        np.where(mean_vai < s.vai_mbl_thr, 1.0 - mean_vai / s.vai_mbl_thr, 0.0) * (1.0 - frac_sno),
        0.0,
    )
    if np.any(((mobile < 0.0) | (mobile > 1.0)) & ~lake):
        raise ValueError("DUSTMod.F90:519 mobile fraction outside [0,1]")
    impossible = (
        np.any(~np.isfinite(gwc))
        or np.any(om_frac < 0.0)
        or np.any((sand < 0.0) | (sand > 100.0))
        or np.any((clay < 0.0) | (clay > 100.0))
        or np.any(np.asarray(h2osoi_vol) < 0.0)
        or np.any(np.asarray(h2osoi_liq) < 0.0)
        or np.any(np.asarray(h2osoi_ice) < 0.0)
        or np.any(np.asarray(rhoa) <= 0.0)
        or np.any(np.asarray(ustar) < 0.0)
        or np.any(np.asarray(forc_pbot) <= 0.0)
        or np.any(np.asarray(forc_t) <= 0.0)
        or np.any(np.asarray(ram1) < 0.0)
        or np.any((np.asarray(land_fraction) < 0.0) | (np.asarray(land_fraction) > 1.0))
        or np.any(weight < 0.0)
    )
    if impossible:
        warnings.warn("CTSM soil input outside physical range", LeungInputWarning)
    flux, diagnostics = dust_emission(
        ustar,
        rhoa,
        gwc,
        mean_vai,
        frac_sno,
        clay,
        rockfrc,
        vegefrc,
        sparfrc,
        roughfct,
        liqfrac=liquid_fraction,
        settings=s,
        mobile_fraction=mobile,
        return_diagnostics=True,
    )
    flux = np.where(active & ~lake, flux, 0.0)
    flux = np.where(np.asarray(land_fraction) < s.coast_cutoff, flux * s.coast_scale, flux)
    setup = init_dust_vars(s) if setup is None else setup
    bins = flux[..., None] * setup["ovr_src_snk_mss"].sum(axis=0)
    dep = dust_dry_dep(
        np.asarray(ustar),
        np.asarray(rhoa),
        np.asarray(forc_pbot),
        np.asarray(forc_t),
        np.asarray(ram1),
        settings=s,
        setup=setup,
    )
    dep["vlc_trb"] = np.where(active[..., None], dep["vlc_trb"], initial_value)
    diagnostics["eta"] = np.where(mobile > 0.0, diagnostics["eta"], 0.0)
    diagnostics["f_eff"] = np.where(
        (mobile > 0.0) & (mean_vai <= s.drag_gate),
        drag_partition(mean_vai, rockfrc, vegefrc, sparfrc, roughfct, settings=s),
        0.0,
    )
    for key in (
        "u_ft0",
        "u_ft",
        "u_it",
        "u_st",
        "Cd",
        "frag_expt",
        "u_s",
        "lnd_frc_mbl",
        "eta",
        "gwc",
    ):
        if key in diagnostics:
            diagnostics[key] = np.where(lake, initial_value, diagnostics[key])
    diagnostics["f_eff"] = np.where(lake, initial_value, diagnostics["f_eff"])
    return {
        "flx_mss_vrt_dst": bins,
        "flx_mss_vrt_dst_tot": np.where(lake, initial_value, bins.sum(axis=-1)),
        "flux_before_split": np.where(lake, initial_value, flux),
        "watsat": watsat,
        "gwc": np.where(lake, initial_value, gwc),
        "liqfrac": np.where(lake, initial_value, liquid_fraction),
        "tlai_lu": mean_vai,
        "gwc_thr": moisture_factor(gwc, clay, settings=s)[1],
        "f_clay": np.minimum(clay / 100.0, s.clay_cap),
        **diagnostics,
        **setup,
        **dep,
    }


# ===========================================================================
# Box-model demo -- run with: python dust_leung2023.py
# ===========================================================================
