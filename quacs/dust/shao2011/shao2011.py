"""
Dust Emission Scheme: UoC Shao 2001/2004/2011 (WRF-Chem dust_opt=4)
======================================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Shao, Y., Ishizuka, M., Mikami, M., and Leys, J. F. (2011).
    Parameterization of size-resolved dust emission and validation with
    measurements. J. Geophys. Res., 116, D08203.
    https://doi.org/10.1029/2010JD014527
    (Shao2011 is the "minimal" / simplified form of Shao 2004; the vertical
     dust flux is Eq. 34 of this paper.)

    Shao, Y. (2004). Simplification of a dust emission scheme and comparison
    with data. J. Geophys. Res., 109, D10202.
    https://doi.org/10.1029/2003JD004372

    Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold
    friction velocity. J. Geophys. Res., 105, 22437-22443.
    (dry-soil threshold u_t0 used here.)

    White, B. R. (1979). Soil transport by winds on Mars. J. Geophys. Res.,
    84, 4643-4651. (streamwise saltation flux Q.)

    Fecan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of
    the increase of the aeolian erosion threshold wind friction velocity due
    to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157.
    (soil-moisture correction to the threshold.)

    Raupach, M. R., Gillette, D. A., and Leys, J. F. (1992). The effect of
    roughness elements on wind erosion threshold. J. Geophys. Res., 98,
    3023-3029. (roughness/drag-partition correction to the threshold.)

    Zeng, Y., Wang, M., Zhao, C., Chen, S., Liu, Z., Huang, X., and Gao, Y.
    (2020). WRF-Chem v3.9 simulations of the East Asian dust storm in May
    2017. Geosci. Model Dev., 13, 2125-2147.
    https://doi.org/10.5194/gmd-13-2125-2020
    (Section 2.2.2 gives the simplified vertical-flux form, Eq. 4, used as the
     scientific summary here; Appendices A and B document the WRF-Chem version
     differences and a post-v3.9 bug fix in this scheme.)

    Fortran reference: WRF-Chem v4.8.0, commit 06d4240
      chem/module_uoc_dust.F   subroutine uoc_dust_driver (erodibility gate,
                               soil-class loop, psd_create call, unit handling)
      chem/module_qf03.F       subroutine qf03_driver (imod=1, 2 and 3 branches,
                               L223-339),
                               ustart0 (Shao-Lu 2000), qwhite (White 1979),
                               h_c (Fecan 1999), r_c (Raupach 1992)
      chem/module_soilpsd.F    minimally/fully-dispersed soil PSD coefficients
                               for 12 soil texture classes
    Scheme authors (Fortran): Y. Shao (U. Cologne); WRF implementation by
    JY Kang, CL Wu, and M. Klose.

QUACS Contract
--------------
Standalone Python model; NumPy is the only runtime dependency. It provides
five-bin rate and full UoC driver APIs; pytest is needed only to run tests.

Pre-API note:
    Once the QUACS scripting API is finalized, the expected signature is:
        @scheme(category='dust_emission')
        def dust_emission_shao2011(state: ModelState, config: Config) -> SurfaceFlux:
    Until then, plain NumPy arrays are used. This script is the scientific
    specification that drives that API design.

How Shao2011 differs from Ginoux (dust_opt=1) and AFWA (dust_opt=3)
-------------------------------------------------------------------
    1. ERODIBILITY IS A LOCATION GATE ONLY. Ginoux multiplies the flux by the
       erodibility field S, so S scales the emission magnitude everywhere.
       Shao2011 permits emission wherever erodibility > 0 and otherwise does
       not use it in the flux at all (see module_uoc_dust.F: the whole cell
       calculation runs only inside `if (sum(erod(i,j,:)).gt.0.)`). This is the
       single most important behavioural difference and is the reason the two
       schemes produce very different spatial patterns (Zeng et al. 2020,
       Sect. 3.2): over the same domain Shao2011 emits more over the Gobi and
       less over the Taklimakan than GOCART.

    2. IT IS DRIVEN BY FRICTION VELOCITY u*, like AFWA, not by 10-m wind.

    3. THE SOIL PARTICLE-SIZE DISTRIBUTION IS THE ENGINE. Each soil texture
       class carries a "minimally dispersed" PSD (aggregated soil as it sits
       in the field) and a "fully dispersed" PSD (after full disaggregation).
       The mass fraction of free dust at size d, eta_mi, IS the
       minimally-dispersed PSD increment dpsdm(d). There is no fixed Kok-style
       emitted-dust distribution as in AFWA; the size spectrum comes from soil
       texture.

    4. SALTATION AND BOMBARDMENT ARE EXPLICIT. The vertical dust flux is the
       streamwise saltation flux Q converted to vertical flux by g/u*^2 and
       amplified by the bombardment factor (1 + sigma_m).

Physics summary
---------------
Per soil texture class c present in a cell (weighted by its area fraction),
the simplified Shao2011 vertical flux for dust size bin j is (Zeng et al. 2020,
Eq. 4; Shao et al. 2011, Eq. 34; qf03 imod=3):

    F(d_j) = (1 - cf) * c_y * eta_mi(d_j) * (1 + sigma_m) * Q * g / u*^2

with:
    c_y       = 1.0e-5                          dimensionless coefficient (cys)
    eta_mi    = dpsdm(d_j)                       minimally-dispersed PSD increment
    sigma_m   = 12 * zeta^2 * (1 + 14 * zeta)    bombardment efficiency
    zeta      = u* * sqrt(rho_s / p)             rho_s = 1000, p = 30000
    cf        = vegetation cover fraction [-]
    Q         = sum_i (1 - cf) * Q_white(u*t_i, u*) * dpsds_i   total saltation
    dpsds     = ghl*dpsdm + (1 - ghl)*dpsdf      sediment PSD
    ghl       = exp(-(u* - u*t_min)^3)            dispersion weight

The dry threshold per sand size (Shao & Lu 2000) is:
    u*t0(d) = sqrt( f * (sigma * g * d + gamma / (rho_a * d)) ),
    f = 0.0123,  gamma = 1.65e-4 N/m,  sigma = rho_p / rho_a
and is corrected for roughness (Raupach 1992) and soil moisture (Fecan 1999):
    u*t = R_rough(lambda) * H_moist(w) * u*t0(d).

The streamwise saltation flux per sand size (White 1979) is:
    Q_white = c * (rho_a / g) * u*^3 * (1 - b) * (1 + b)^2,  b = u*t/u*,
    c = 2.3   (note: the qf03 comment cites 2.6, but the active code uses 2.3.)

Emission is zero where erodibility = 0, and zero for a soil class when
u* <= u*t_min for that class.

Inputs (dynamic, per time step from MPAS-A):
    ustar : np.ndarray [x, y]     friction velocity u* [m/s]
    gwet  : np.ndarray [x, y]     volumetric soil moisture [-], in [0, 1]
    cf    : np.ndarray [x, y]     vegetation cover [%]

Inputs (static, loaded once at initialization):
    erod      : np.ndarray [x, y, 3]     erodibility (sum > 0 is a source)
    soil_frac : np.ndarray [x, y, 12]    area fraction of each soil texture class
    rhoa      : np.ndarray [x, y]        near-surface air density [kg/m^3]
    The per-class minimally/fully-dispersed PSDs are built once by initialize()
    from the soil-texture coefficient tables below (no external files needed).

Configurable:
    c_y   : float   emission coefficient cys (default 1.0e-5, qf03 L130).
    sf_surface_physics: land model dry-soil table (default 2, Noah).
    imod: 1 (Shao 2001), 2 (Shao 2004), or 3 (Shao 2011; default).
    The White saltation constant c (default 2.3) and the soil/particle
    constants are module-level and source-pinned.

Output:
    F_p : np.ndarray [x, y, n_bins]   surface emission flux [kg m^-2 s^-1]
                                      over the 5 GOCART dust bins (upper-cut
                                      diameters 2, 3.6, 6, 12, 20 micron).
                                      QUACS runtime converts to tendency:
                                      dx/dt = F_p / (rho_air * dz)
    uoc_dust_driver also returns updated ug/kg tracers, accumulated kg/m2
    emissions, top-layer moisture and three WRF diagnostics.

Notes on faithfulness and source quirks (flag for methods documentation)
------------------------------------------------------------------------
  * The (1 - cf) factor is applied TWICE in the WRF implementation: once inside
    the saltation flux q (module_uoc_dust path / qf03 q = (1-cf)*qwhite) and
    again on the final dust flux fff = (1-cf)*fff. This script reproduces the
    source exactly, so the net cover dependence is (1 - cf)^2. This holds for
    imod=1,2,3 in qf03, so it appears intentional in the implementation rather
    than specific to Shao2011, but it is worth confirming with the team.
  * White constant: code value 2.3 is used; the inline comment in qwhite cites
    2.6. Pinned to the active code value.
  * The size grid and soil PSD coefficients follow the pinned Fortran files;
    NumPy vectorizes psd_create and uses float64 throughout.
  * Zeng et al. (2020) corrected a post-v3.9 bug in this scheme; this script is
    pinned to WRF v4.8.0 (06d4240). Exact reproduction of the paper's numbers
    requires the v3.9-plus-fix variant described in their Appendices A and B.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, fields

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants and scheme parameters
# Source: chem/module_qf03.F and chem/module_uoc_dust.F (WRF v4.8.0, 06d4240)
# ---------------------------------------------------------------------------
GRAV = 9.81  # m/s^2


class ShaoDustWarning(UserWarning):
    """Warnings for active UoC dust cells with unusual inputs or diagnostics."""


# Emission coefficient (qf03.F L130: cys = 0.00001)
CYS_DEFAULT = 1.0e-5

# Soil / particle constants (qf03_driver locals)
RHOS = 1000.0  # soil bulk density rho_s [kg/m^3]   (qf03.F: rhos = 1000)
PHL = 30000.0  # plastic pressure p [N/m^2]         (qf03.F: phl = 30000)
RHOP_DEFAULT = 2650.0  # particle density rho_p [kg/m^3]    (uoc_dust.F: rhop = 2650)

# Shao & Lu (2000) threshold coefficients (qf03.F: ustart0)
F_SHAO = 0.0123  # dimensionless
GAMMA_SHAO = 1.65e-4  # N/m

# White (1979) saltation constant (qf03.F: qwhite, c = 2.3; comment cites 2.6)
C_WHITE = 2.3

# Raupach (1992) roughness correction (qf03.F: r_c)
C_LAMBDA = 0.35  # frontal-area-index coefficient (qf03_driver: c_lambda)
RAUP_SIG = 1.0  # basal-to-frontal area ratio
RAUP_M = 0.5  # parameter < 1
RAUP_BETA = 200.0  # drag-coefficient ratio (Shao recommendation)
COVER_CAP = 0.95
GHL_POWER = 3.0
SIGMA_SCALE, SIGMA_SLOPE = 12.0, 14.0
IMPACT_ANGLE, IMPACT_WIND, IMPACT_CUBIC = 13.0, 10.0, 7.5

# Size grid (psd_create in module_uoc_dust.F)
IMAX = 100  # number of size intervals
JMAX = 4  # log-normal modes per soil PSD
DREF = 1000.0  # reference particle size [micron]
PHI_UP = 10.0  # phi upper
PHI_LO = -1.0  # phi lower
DCUT = 20.0  # dust cutoff diameter [micron]

# GOCART output bins: upper-cut diameters [micron] (qf03_driver: dbin)
DBIN = np.array([2.0, 3.6, 6.0, 12.0, 20.0])
N_BINS = len(DBIN)

# ---------------------------------------------------------------------------
# Soil texture classes (1..12) in WRF/Shao order
#   1 sand, 2 loamy sand, 3 sandy loam, 4 silt loam, 5 silt, 6 loam,
#   7 sandy clay loam, 8 silty clay loam, 9 clay loam, 10 sandy clay,
#   11 silty clay, 12 clay
# ---------------------------------------------------------------------------
SOIL_CLASS_NAMES = [
    "sand",
    "loamy_sand",
    "sandy_loam",
    "silt_loam",
    "silt",
    "loam",
    "sandy_clay_loam",
    "silty_clay_loam",
    "clay_loam",
    "sandy_clay",
    "silty_clay",
    "clay",
]
N_SOIL = len(SOIL_CLASS_NAMES)

# Residual / dry soil moisture per soil class (uoc_dust.F: data thr/.../)
# Used as wr in the Fecan correction. Indexed by WRF soil class 1..12.
WR_THR = np.array(
    [0.001, 0.003, 0.037, 0.061, 0.072, 0.049, 0.084, 0.110, 0.095, 0.126, 0.141, 0.156]
)
NOAH_DRYSMC = np.array(
    [0.010, 0.028, 0.047, 0.084, 0.061, 0.066, 0.069, 0.120, 0.103, 0.100, 0.126, 0.138]
)
RUC_DRYSMC = np.array(
    [0.002, 0.035, 0.041, 0.034, 0.034, 0.050, 0.068, 0.060, 0.050, 0.070, 0.070, 0.068]
)

# Fecan (1999) coefficients (qf03.F: h_c), indexed by Shao class isl 1..12
FECAN_A = np.array(
    [21.19, 33.03, 44.87, 17.79, 20.81, 23.83, 26.84, 29.86, 27.51, 25.17, 22.82, 20.47]
)
FECAN_B = np.array([0.68, 0.71, 0.85, 0.61, 0.66, 0.71, 0.75, 0.80, 0.75, 0.70, 0.64, 0.59])

# WRF soil class (cc, 1..12) -> Shao class (isl, 1..12) for Fecan a,b
# (uoc_dust.F mapping). Stored 0-indexed: CC2ISL[cc-1] = isl-1.
_CC2ISL_1BASED = {1: 1, 2: 2, 3: 3, 4: 5, 5: 6, 6: 4, 7: 7, 8: 9, 9: 8, 10: 10, 11: 11, 12: 12}
CC2ISL = np.array([_CC2ISL_1BASED[cc] - 1 for cc in range(1, N_SOIL + 1)])

# ---------------------------------------------------------------------------
# Soil PSD log-normal coefficients (module_soilpsd.F)
# Each entry is 4 modes of (w, dln, sig): weight, mean log-diameter, log-sigma.
# "_M" = minimally dispersed, "_F" = fully dispersed.
# Transcribed verbatim (Fortran data fills (3,4) column-major, so each triple
# below is one mode).
# ---------------------------------------------------------------------------
SOIL_PSD_M = {
    "sand": [
        (0.0287, 3.6153, 0.2775),
        (0.2811, 4.9918, 0.3023),
        (0.0516, 3.9315, 0.1417),
        (0.6387, 4.7173, 0.2432),
    ],
    "loamy_sand": [
        (0.5000, 6.0674, 0.4039),
        (0.1997, 4.3282, 0.3998),
        (0.2191, 5.2793, 0.3488),
        (0.0812, 6.6426, 0.2216),
    ],
    "sandy_loam": [
        (0.0240, 6.6566, 0.1000),
        (0.0536, 6.0663, 0.1227),
        (0.3184, 5.1840, 0.7462),
        (0.6039, 6.0685, 0.4063),
    ],
    "silt_loam": [
        (0.0278, 5.2068, 0.1921),
        (0.5000, 4.3275, 0.4544),
        (0.3054, 3.5848, 1.0721),
        (0.1669, 4.1432, 0.1877),
    ],
    "silt": [
        (0.0278, 5.2068, 0.1921),
        (0.5000, 4.3275, 0.4544),
        (0.3054, 3.5848, 1.0721),
        (0.1669, 4.1432, 0.1877),
    ],
    "loam": [
        (0.0695, 6.9010, 0.1000),
        (0.1047, 6.6666, 0.1614),
        (0.5000, 5.5720, 0.7295),
        (0.3258, 6.3005, 0.3549),
    ],
    "sandy_clay_loam": [
        (0.4999, 5.1720, 0.3064),
        (0.2490, 4.6158, 0.2783),
        (0.0139, 4.9110, 0.1000),
        (0.2372, 5.0185, 0.9259),
    ],
    "silty_clay_loam": [
        (1.2597, 4.7986, 0.3751),
        (0.8107, 5.2549, 0.3047),
        (0.4482, 5.1246, 1.2550),
        (0.0000, 0.0000, 0.0000),
    ],
    "clay_loam": [
        (0.1842, 6.3110, 0.2071),
        (0.4243, 6.0792, 0.4049),
        (0.3273, 5.5946, 0.7726),
        (0.0642, 6.5793, 0.1000),
    ],
    "sandy_clay": [
        (0.3124, 4.1426, 0.1717),
        (0.9564, 3.9501, 1.7750),
        (1.0340, 4.3088, 0.4340),
        (0.0000, 0.0000, 0.0000),
    ],
    "silty_clay": [
        (0.2709, 4.9160, 0.1971),
        (0.0631, 4.5807, 0.1635),
        (0.1350, 3.8960, 0.8092),
        (0.5310, 4.5301, 0.4887),
    ],
    "clay": [
        (0.2408, 4.5855, 0.6331),
        (0.0594, 3.3126, 1.1665),
        (0.0273, 5.3894, 0.1000),
        (0.6725, 5.3148, 0.3924),
    ],
}
SOIL_PSD_F = {
    "sand": [
        (0.0231, 3.6724, 0.2341),
        (0.0362, 3.9598, 0.1257),
        (0.2628, 4.9933, 0.2986),
        (0.6779, 4.7374, 0.2498),
    ],
    "loamy_sand": [
        (0.1354, 5.5976, 0.4288),
        (0.1073, 2.3499, 1.0898),
        (0.1692, 4.0550, 0.2113),
        (0.5880, 4.1982, 0.7748),
    ],
    "sandy_loam": [
        (0.0115, 5.3900, 0.1141),
        (0.3043, 4.6980, 1.0132),
        (0.0840, 6.6115, 0.1339),
        (0.6002, 6.0494, 0.3178),
    ],
    "silt_loam": [
        (0.1816, 3.1175, 1.0169),
        (0.4454, 4.3491, 0.6154),
        (0.0568, 4.1250, 0.1000),
        (0.3162, 4.1594, 0.3017),
    ],
    "silt": [
        (0.1816, 3.1175, 1.0169),
        (0.4454, 4.3491, 0.6154),
        (0.0568, 4.1250, 0.1000),
        (0.3162, 4.1594, 0.3017),
    ],
    "loam": [
        (0.0378, 5.0205, 0.5601),
        (0.0511, 0.5580, 0.3886),
        (0.4003, 1.7677, 0.6877),
        (0.5108, 2.9973, 0.5489),
    ],
    "sandy_clay_loam": [
        (0.1364, 3.3869, 1.3277),
        (0.0642, 6.1715, 0.3463),
        (0.0767, 4.4147, 0.2243),
        (0.7227, 4.9219, 0.4983),
    ],
    "silty_clay_loam": [
        (0.5844, 4.6079, 0.6141),
        (0.3304, 5.2050, 0.2897),
        (0.0522, 7.0553, 1.0000),
        (0.0330, 0.6931, 1.0000),
    ],
    "clay_loam": [
        (0.3988, 2.0568, 0.9444),
        (0.5000, 3.5758, 0.6679),
        (0.0145, 4.4626, 0.1000),
        (0.0867, 5.0406, 0.4590),
    ],
    "sandy_clay": [
        (1.1285, 4.2288, 0.4296),
        (0.7275, 3.8558, 1.2898),
        (0.2438, 4.1222, 0.1379),
        (0.0000, 0.0000, 0.0000),
    ],
    "silty_clay": [
        (0.3927, 3.6159, 0.5702),
        (0.1160, 1.3131, 0.2704),
        (0.4588, 2.3768, 0.6913),
        (0.0325, 1.0150, 0.1154),
    ],
    "clay": [
        (0.5000, 3.0166, 0.7156),
        (0.1124, 0.6297, 0.4159),
        (0.3698, 1.5791, 0.6059),
        (0.0178, 3.8284, 0.1509),
    ],
}


@dataclass(frozen=True)
class Shao2011Settings:
    c_y: float = field(
        default=CYS_DEFAULT, metadata={"source": "qf03.F:130", "what": "vertical dust scale"}
    )
    rhos: float = field(default=RHOS, metadata={"source": "qf03.F:128", "what": "soil density"})
    phl: float = field(default=PHL, metadata={"source": "qf03.F:129", "what": "plastic pressure"})
    rhop: float = field(
        default=RHOP_DEFAULT, metadata={"source": "uoc:81", "what": "particle density"}
    )
    f_shao: float = field(
        default=F_SHAO, metadata={"source": "qf03:382", "what": "dry threshold scale"}
    )
    gamma_shao: float = field(
        default=GAMMA_SHAO, metadata={"source": "qf03:381", "what": "dry threshold cohesion"}
    )
    c_white: float = field(
        default=C_WHITE, metadata={"source": "qf03:408", "what": "saltation scale"}
    )
    c_lambda: float = field(
        default=C_LAMBDA, metadata={"source": "qf03:88", "what": "frontal area coefficient"}
    )
    raup_sig: float = field(
        default=RAUP_SIG, metadata={"source": "qf03:521", "what": "drag basal area factor"}
    )
    raup_m: float = field(
        default=RAUP_M, metadata={"source": "qf03:521", "what": "drag partition factor"}
    )
    raup_beta: float = field(
        default=RAUP_BETA, metadata={"source": "qf03:521", "what": "drag coefficient ratio"}
    )
    cover_cap: float = field(
        default=COVER_CAP, metadata={"source": "qf03:138", "what": "vegetation cover cap"}
    )
    imax: int = field(default=IMAX, metadata={"source": "uoc:71", "what": "size grid intervals"})
    jmax: int = field(default=JMAX, metadata={"source": "uoc:72", "what": "PSD mode count"})
    dref: float = field(default=DREF, metadata={"source": "uoc:306", "what": "size grid reference"})
    phi_up: float = field(
        default=PHI_UP, metadata={"source": "uoc:323", "what": "size grid upper phi"}
    )
    phi_lo: float = field(
        default=PHI_LO, metadata={"source": "uoc:323", "what": "size grid lower phi"}
    )
    dcut: float = field(default=DCUT, metadata={"source": "uoc:80", "what": "dust size cutoff"})
    dbin: tuple = field(
        default_factory=lambda: tuple(DBIN),
        metadata={"source": "qf03:102", "what": "dust bin cutoffs"},
    )
    wr_thr: tuple = field(
        default_factory=lambda: tuple(WR_THR),
        metadata={"source": "uoc:99", "what": "WRF residual moisture"},
    )
    noah_drysmc: tuple = field(
        default_factory=lambda: tuple(NOAH_DRYSMC),
        metadata={
            "source": "SOILPARM.TBL:4-15 (STAS table, DRYSMC column, read by noahlsm.F)",
            "what": "Noah residual moisture",
        },
    )
    ruc_drysmc: tuple = field(
        default_factory=lambda: tuple(RUC_DRYSMC),
        metadata={
            "source": "SOILPARM.TBL:26-37 (STAS-RUC table, DRYSMC column, read by ruclsm.F:7324)",
            "what": "RUC residual moisture",
        },
    )
    fecan_a: tuple = field(
        default_factory=lambda: tuple(FECAN_A),
        metadata={"source": "qf03:469", "what": "Fecan moisture scale"},
    )
    fecan_b: tuple = field(
        default_factory=lambda: tuple(FECAN_B),
        metadata={"source": "qf03:470", "what": "Fecan moisture exponent"},
    )
    soil_psd_m: tuple = field(
        default_factory=lambda: tuple(
            tuple(map(tuple, SOIL_PSD_M[name])) for name in SOIL_CLASS_NAMES
        ),
        metadata={"source": "module_soilpsd.F:7-77", "what": "minimal soil size distributions"},
    )
    soil_psd_f: tuple = field(
        default_factory=lambda: tuple(
            tuple(map(tuple, SOIL_PSD_F[name])) for name in SOIL_CLASS_NAMES
        ),
        metadata={"source": "module_soilpsd.F:79-150", "what": "full soil size distributions"},
    )
    ghl_power: float = field(
        default=GHL_POWER, metadata={"source": "qf03:197", "what": "dispersion mix exponent"}
    )
    sigma_scale: float = field(
        default=SIGMA_SCALE, metadata={"source": "qf03:274", "what": "bombardment scale"}
    )
    sigma_slope: float = field(
        default=SIGMA_SLOPE, metadata={"source": "qf03:319", "what": "bombardment slope"}
    )
    impact_angle: float = field(
        default=IMPACT_ANGLE, metadata={"source": "qf03:229", "what": "impact angle degrees"}
    )
    impact_wind: float = field(
        default=IMPACT_WIND, metadata={"source": "qf03:228", "what": "impact wind scale"}
    )
    impact_cubic: float = field(
        default=IMPACT_CUBIC, metadata={"source": "qf03:440", "what": "impact cubic scale"}
    )

    def __post_init__(self):
        for name in (
            "dbin",
            "wr_thr",
            "noah_drysmc",
            "ruc_drysmc",
            "fecan_a",
            "fecan_b",
            "soil_psd_m",
            "soil_psd_f",
        ):
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                value = value.copy()
                value.setflags(write=False)
                object.__setattr__(self, name, value)


def tunable_names(settings_class):
    return tuple(item.name for item in fields(settings_class))


# ---------------------------------------------------------------------------
# Soil PSD construction: psd_create (module_uoc_dust.F)
# ---------------------------------------------------------------------------
def psd_create(coeffs, settings=None):
    """
    Build a discretized particle-size distribution from log-normal modes.

    Faithful port of psd_create (module_uoc_dust.F). Returns the size-interval
    boundaries, midpoints, and the normalized PSD on the standard phi grid.

    Parameters
    ----------
    coeffs : list of (w, dln, sig)   up to JMAX log-normal modes.

    Returns
    -------
    d   : np.ndarray [IMAX+1]   interval boundaries [micron]
    dm  : np.ndarray [IMAX]     interval midpoints (delta-d / delta-ln-d) [micron]
    psd : np.ndarray [IMAX]     density p(d) [micron^-1]
    dpsd: np.ndarray [IMAX]     mass fraction increment Delta P per interval [-]
    ppsd: np.ndarray [IMAX]     cumulative P(d) [-]
    """
    settings = Shao2011Settings() if settings is None else settings
    cn = 1.0 / np.sqrt(2.0 * 3.14159)
    i = np.arange(0, settings.imax + 1)
    phi = settings.phi_up - i * (settings.phi_up - settings.phi_lo) / settings.imax
    d = settings.dref / 2.0**phi  # [IMAX+1], increasing
    dm = (d[1:] - d[:-1]) / np.log(d[1:] / d[:-1])  # [IMAX]

    pp = np.zeros(settings.imax)
    for w, dln, sig in coeffs[: settings.jmax]:
        if (w > 1.0e-7) and (sig != 0.0):
            pp += w * cn / sig * np.exp(-((np.log(dm) - dln) ** 2) / (2.0 * sig**2))

    dpsd = pp * (np.log(d[1:]) - np.log(d[:-1]))
    ppsd = np.cumsum(dpsd)
    psd = pp / dm

    norm = ppsd[-1]
    dpsd = dpsd / norm
    psd = psd / norm
    ppsd = ppsd / norm
    return d, dm, psd, dpsd, ppsd


# ---------------------------------------------------------------------------
# Threshold friction velocity: Shao & Lu (2000) (qf03.F: ustart0)
# ---------------------------------------------------------------------------
def threshold_velocity_shao_lu(d_micron, sigma, rhoa, settings=None):
    """
    Dry, bare-surface threshold friction velocity per particle size.

    u*t0(d) = sqrt( f * (sigma * g * d + gamma / (rhoa * d)) )

    Parameters
    ----------
    d_micron : np.ndarray [n_size]   particle diameter [micron]
    sigma    : float                 particle-to-air density ratio rho_p/rho_a
    rhoa     : float                 air density [kg/m^3]

    Returns
    -------
    u_t0 : np.ndarray [n_size]   threshold friction velocity [m/s]
    """
    settings = Shao2011Settings() if settings is None else settings
    dm = d_micron * 1.0e-6  # micron -> m
    val = settings.f_shao * (sigma * GRAV * dm + settings.gamma_shao / (rhoa * dm))
    return np.sqrt(val)


# ---------------------------------------------------------------------------
# Streamwise saltation flux: White (1979) (qf03.F: qwhite)
# ---------------------------------------------------------------------------
def saltation_flux_white(u_t, ustar, rhoa, settings=None):
    """
    White (1979) streamwise sand flux per size, zero below threshold.

    Q = c * (rhoa / g) * u*^3 * (1 - b) * (1 + b)^2,  b = u_t / u*

    u_t and the broadcast result follow array shapes; emission is zero where
    ustar < u_t.
    """
    settings = Shao2011Settings() if settings is None else settings
    b = u_t / ustar
    q = settings.c_white * (rhoa / GRAV) * ustar**3 * (1.0 - b) * (1.0 + b) ** 2
    return np.where(ustar < u_t, 0.0, q)


# ---------------------------------------------------------------------------
# Roughness correction: Raupach et al. (1992) (qf03.F: r_c)
# ---------------------------------------------------------------------------
def roughness_correction_raupach(cf, settings=None):
    """
    Roughness/drag-partition correction R = u*tR / u*tS as a function of
    vegetation cover fraction cf (via frontal area index lambda).

        lambda = -C_LAMBDA * ln(1 - cf)
        R      = sqrt(1 - sig*m*lambda) * sqrt(1 + m*beta*lambda)

    cf is clipped at 0.95 (as in the Fortran) to keep lambda finite. Fully
    covered cells (lambda >= 1/(sig*m)) return a large R (effectively no
    emission), matching r_c returning 999.
    """
    settings = Shao2011Settings() if settings is None else settings
    cf = np.minimum(cf, settings.cover_cap)
    lam = -settings.c_lambda * np.log(1.0 - cf)
    xc = 1.0 / (settings.raup_sig * settings.raup_m)
    inside = 1.0 - settings.raup_sig * settings.raup_m * lam
    r = np.where(
        lam >= xc,
        999.0,
        np.sqrt(np.maximum(inside, 0.0))
        * np.sqrt(1.0 + settings.raup_m * settings.raup_beta * lam),
    )
    return r


# ---------------------------------------------------------------------------
# Soil-moisture correction: Fecan et al. (1999) (qf03.F: h_c)
# ---------------------------------------------------------------------------
def moisture_correction_fecan(w, wr, isl, settings=None):
    """
    Fecan (1999) threshold increase due to soil moisture.

        H = 1                              if w <= wr
        H = sqrt(1 + a(isl) * (w - wr)^b(isl))   otherwise

    Parameters
    ----------
    w   : np.ndarray            volumetric soil moisture [-]
    wr  : float                 residual moisture for the soil class [-]
    isl : int                   Shao soil class index (0-based) for a, b
    """
    settings = Shao2011Settings() if settings is None else settings
    a = settings.fecan_a[isl]
    b = settings.fecan_b[isl]
    excess = np.maximum(w - wr, 0.0)
    return np.where(w <= wr, 1.0, np.sqrt(1.0 + a * excess**b))


# ---------------------------------------------------------------------------
# Initialization: build the static soil PSDs and bin index map once
# ---------------------------------------------------------------------------
def initialize(rhoa_ref=1.225, settings=None):
    """
    Build the static per-soil-class particle-size distributions and the
    GOCART bin index map. Called once before the time loop.

    No external files are required: the soil PSD coefficients are embedded
    (module_soilpsd.F). rhoa_ref is validated for the reporting interface;
    the timestep function uses the per-cell air density.

    Returns
    -------
    psd : dict with arrays shared across all soil classes and the per-class
          dispersed distributions:
            d       [IMAX+1]            interval boundaries [micron]
            dm      [IMAX]              interval midpoints [micron]
            dpsd_m  [N_SOIL, IMAX]      minimally-dispersed increments
            dpsd_f  [N_SOIL, IMAX]      fully-dispersed increments
            idst    int                 number of dust-size intervals (dm <= DCUT)
            bin_hi  [N_BINS]            inclusive upper interval index per bin
    """
    settings = Shao2011Settings() if settings is None else settings
    d, dm, _, _, _ = psd_create(settings.soil_psd_m[0], settings)

    psd_m = np.zeros((N_SOIL, settings.imax))
    psd_f = np.zeros((N_SOIL, settings.imax))
    dpsd_m = np.zeros((N_SOIL, settings.imax))
    dpsd_f = np.zeros((N_SOIL, settings.imax))
    for k, name in enumerate(SOIL_CLASS_NAMES):
        _, _, pm, dpm, _ = psd_create(settings.soil_psd_m[k], settings)
        _, _, pf, dpf, _ = psd_create(settings.soil_psd_f[k], settings)
        psd_m[k] = pm
        psd_f[k] = pf
        dpsd_m[k] = dpm
        dpsd_f[k] = dpf

    # number of dust-size intervals (dm <= DCUT), as in uoc_dust.F idst
    idst = int(np.count_nonzero(dm <= settings.dcut))

    # bin upper-index map: smallest interval index i (1-based on dm) whose
    # boundary d0[i] >= DBIN[n]; reproduces qf03_driver ibin via boundaries d.
    bin_hi = np.zeros(N_BINS, dtype=int)
    for n in range(N_BINS):
        idx = np.where(d[1:] >= settings.dbin[n])[0]  # boundaries d[1..IMAX]
        if not idx.size:
            raise ValueError("[dust_shao2011] wrong dust classes")
        bin_hi[n] = idx[0]

    if not (0.0 <= rhoa_ref <= 5.0):
        raise ValueError(f"[dust_shao2011] rhoa_ref out of range: {rhoa_ref}")

    return dict(
        d=d, dm=dm, psd_m=psd_m, psd_f=psd_f, dpsd_m=dpsd_m, dpsd_f=dpsd_f, idst=idst, bin_hi=bin_hi
    )


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def _drysmc(sf_surface_physics, active, settings):
    if sf_surface_physics == 3:
        return settings.ruc_drysmc
    if sf_surface_physics in (0, 1, 5, 7, 8):
        return settings.wr_thr
    if sf_surface_physics == 4:
        # WRF commit 49e5342870e9fb3099a247d3cfc04917a7441bc4.
        return settings.noah_drysmc
    if sf_surface_physics not in (-1, 2) and active.any():
        warnings.warn(
            "[dust_shao2011] unknown sf_surface_physics; using Noah DRYSMC", ShaoDustWarning
        )
    return settings.noah_drysmc


def _flux_and_diagnostics(
    ustar,
    gwet,
    cf,
    erod,
    soil_frac,
    rhoa,
    psd,
    sf_surface_physics,
    imod,
    isltyp=None,
    settings=None,
):
    """Evaluate Q:115-359 once and return rate, capped moisture and diagnostics."""
    settings = Shao2011Settings() if settings is None else settings
    ustar, gwet, cf = (np.asarray(value, dtype=float) for value in (ustar, gwet, cf))
    erod, soil_frac, rhoa = (np.asarray(value, dtype=float) for value in (erod, soil_frac, rhoa))
    shape = ustar.shape
    if ustar.ndim != 2 or any(np.shape(a) != shape for a in (gwet, cf, rhoa)):
        raise ValueError("[dust_shao2011] inconsistent surface grid shapes")
    if erod.shape != (*shape, 3):
        raise ValueError("[dust_shao2011] erod must have three parts on the last axis")
    if soil_frac.shape != (*shape, N_SOIL):
        raise ValueError("[dust_shao2011] soil_frac must have 12 classes")
    if imod not in (1, 2, 3):
        raise ValueError("[dust_shao2011] imod must be 1, 2 or 3")
    if isltyp is not None and np.shape(isltyp) != shape:
        raise ValueError("[dust_shao2011] isltyp has wrong shape")
    active = erod.sum(axis=-1) > 0.0
    if np.any((erod[active] < 0.0) | ~np.isfinite(erod[active])):
        warnings.warn(
            "[dust_shao2011] active erod has negative or nonfinite parts", ShaoDustWarning
        )
    if np.any((cf[active] < 0.0) | (cf[active] > 100.0) | ~np.isfinite(cf[active])):
        warnings.warn(
            "[dust_shao2011] active cover outside 0..100 percent or nonfinite", ShaoDustWarning
        )
    if np.any((soil_frac[active] < 0.0) | ~np.isfinite(soil_frac[active])):
        warnings.warn(
            "[dust_shao2011] active soil fraction is negative or nonfinite", ShaoDustWarning
        )
    total_soil = np.zeros(shape)
    for c in range(N_SOIL):
        total_soil = total_soil + soil_frac[:, :, c]
    if np.any(active & ((total_soil <= 0.0) | ~np.isfinite(total_soil))):
        warnings.warn("[dust_shao2011] active soil fraction sum is not positive", ShaoDustWarning)
    if np.any(active & (soil_frac != 0.0).any(axis=-1) & (gwet < 0.0)):
        raise ValueError("[dust_shao2011] soil moisture correction (h_c): w < 0")
    if np.any(active & (gwet > 1.0)):
        warnings.warn("[dust_shao2011] active moisture > 1 reset to 1", ShaoDustWarning)
    if np.any(active & ~np.isfinite(gwet)):
        warnings.warn("[dust_shao2011] active moisture is nonfinite", ShaoDustWarning)
    changed_wet = np.where(active, np.minimum(gwet, 1.0), gwet)
    if np.any(active & ((rhoa <= 0.0) | (ustar < 0.0) | ~np.isfinite(rhoa) | ~np.isfinite(ustar))):
        warnings.warn(
            "[dust_shao2011] active density or friction velocity is impossible", ShaoDustWarning
        )
    wr = _drysmc(sf_surface_physics, active, settings)
    cf_fraction = np.where(active, cf / 100.0, 0.0)
    cf_fraction = np.minimum(cf_fraction, settings.cover_cap)
    with np.errstate(all="ignore"):
        rough = roughness_correction_raupach(cf_fraction, settings)
    one_minus_cf = 1.0 - cf_fraction
    wind = np.where(ustar > 0.0, ustar, 1.0e-12)
    density = np.where(active, rhoa, 1.0)
    dm = psd["dm"]
    idst = psd["idst"]
    bin_hi = psd["bin_hi"]
    dm_m = dm * 1.0e-6
    with np.errstate(all="ignore"):
        u_t0 = np.sqrt(
            settings.f_shao
            * (
                (settings.rhop / density)[:, :, None] * GRAV * dm_m
                + settings.gamma_shao / (density[:, :, None] * dm_m)
            )
        )
    fff = np.zeros((*shape, idst))
    ust_diag = np.full(shape, np.nan)
    rough_diag = np.full(shape, np.nan)
    moist_diag = np.full(shape, np.nan)
    for c in range(N_SOIL):
        fraction = soil_frac[:, :, c]
        present = active & (fraction != 0.0)
        if not present.any():
            continue
        hmoist = moisture_correction_fecan(changed_wet, wr[c], CC2ISL[c], settings)
        u_t = rough[:, :, None] * hmoist[:, :, None] * u_t0
        u_t_min = np.minimum(999.0, u_t.min(axis=-1))
        ust_diag[present] = u_t_min[present]
        rough_diag[present] = rough[present]
        if isltyp is not None:
            dominant = present & (np.asarray(isltyp) == c + 1)
            moist_diag[dominant] = hmoist[dominant]
        emit = present & ~(ustar <= u_t_min)
        if not emit.any():
            # Q:349 applies the class weight even when Q:188-195 set fbin=0.
            with np.errstate(all="ignore"):
                weight = fraction / total_soil
                fff += np.where(present[:, :, None], weight[:, :, None] * 0.0, 0.0)
            continue
        with np.errstate(all="ignore"):
            q = one_minus_cf[:, :, None] * saltation_flux_white(
                u_t, wind[:, :, None], density[:, :, None], settings
            )
            ghl = np.exp(-(np.where(emit, ustar - u_t_min, 0.0) ** settings.ghl_power))
            dpsdm = psd["dpsd_m"][c]
            dpsdf = psd["dpsd_f"][c]
            dpsds = ghl[:, :, None] * dpsdm + (1.0 - ghl[:, :, None]) * dpsdf
            ffq = q * dpsds
            if imod == 3:
                zeta = wind * np.sqrt(settings.rhos / settings.phl)
                sigma_m = settings.sigma_scale * zeta**2 * (1.0 + settings.sigma_slope * zeta)
                a1 = settings.c_y * dpsdm[:idst][None, None, :]
                a3 = (ffq.sum(axis=-1) * GRAV / wind**2)[:, :, None]
                rate = one_minus_cf[:, :, None] * a1 * (1.0 + sigma_m[:, :, None]) * a3
            elif imod == 2:
                a3 = GRAV / wind**2
                zeta = wind * np.sqrt(settings.rhos / settings.phl)
                sigma_m = settings.sigma_scale * zeta**2 * (1.0 + settings.sigma_slope * zeta)
                ratio = np.minimum(psd["psd_m"][c, :idst] / psd["psd_f"][c, :idst], 1.0)
                a1 = (
                    settings.c_y
                    * dpsdf[:idst]
                    * ((1.0 - ghl[:, :, None]) + ghl[:, :, None] * ratio)
                )
                rate = (
                    one_minus_cf[:, :, None]
                    * a1
                    * (1.0 + sigma_m[:, :, None])
                    * (ffq[:, :, idst:].sum(axis=-1) * a3)[:, :, None]
                )
            else:
                a3 = GRAV / wind**2
                ratio = np.minimum(psd["psd_m"][c, :idst] / psd["psd_f"][c, :idst], 1.0)
                a1 = settings.c_y * ((1.0 - ghl[:, :, None]) + ghl[:, :, None] * ratio)
                rate = np.zeros((*shape, idst))
                alpha = settings.impact_angle * 3.14159 / 180.0
                u0 = settings.impact_wind * wind
                for i in range(idst, settings.imax):
                    diameter = dm_m[i]
                    mass = 3.1415927 * settings.rhop * diameter**3 / 6.0
                    beta = np.sqrt(settings.phl * 2.0 * diameter / mass)
                    t1 = u0**2 / beta**2 * (np.sin(2.0 * alpha) - 4.0 * np.sin(alpha) ** 2)
                    t2 = u0 * np.sin(alpha) / beta
                    t3 = settings.impact_cubic * 3.1415927 * t2**3 / diameter
                    omega = diameter * (t1 + t3)
                    a2 = ffq[:, :, i] * a3 / mass
                    a_dust = dpsdf[:idst] * settings.rhos * omega[:, :, None]
                    a_dust = a_dust + np.where(
                        dpsdf[:idst] >= dpsdm[:idst], (dpsdf[:idst] - dpsdm[:idst]) * mass, 0.0
                    )
                    rate += a1 * a2[:, :, None] * a_dust
                rate *= one_minus_cf[:, :, None]
            weight = fraction / total_soil
            fff += np.where(
                present[:, :, None], weight[:, :, None] * np.where(emit[:, :, None], rate, 0.0), 0.0
            )
    if np.any(active & (soil_frac == 0.0).all(axis=-1)):
        warnings.warn(
            "[dust_shao2011] all 12 soil fractions are zero; diagnostics undefined", ShaoDustWarning
        )
    if isltyp is not None and np.any(
        active & np.isnan(moist_diag) & (soil_frac != 0.0).any(axis=-1)
    ):
        warnings.warn(
            "[dust_shao2011] smois_cor undefined for zero-fraction isltyp", ShaoDustWarning
        )
    flux = np.zeros((*shape, N_BINS))
    lo = 0
    for n, hi in enumerate(bin_hi):
        flux[:, :, n] = fff[:, :, lo : hi + 1].sum(axis=-1)
        lo = hi + 1
    ust_diag[~active] = -9999.0
    rough_diag[~active] = 1.0
    moist_diag[~active] = 1.0
    return flux, changed_wet, ust_diag, rough_diag, moist_diag


def dust_emission_shao2011(
    ustar, gwet, cf, erod, soil_frac, rhoa, psd, sf_surface_physics=2, imod=3, settings=None
):
    """Return five-bin UoC dust rates [kg m-2 s-1] for WRF cover percent.

    `erod` has three parts on the last axis. `bems` in WRF is this rate times
    `dt` [kg m-2]. `imod` selects Shao 2001, 2004 or 2011 (1, 2 or 3).
    """
    return _flux_and_diagnostics(
        ustar, gwet, cf, erod, soil_frac, rhoa, psd, sf_surface_physics, imod, settings=settings
    )[0]


def uoc_dust_driver(
    dust,
    rho_phy,
    dz8w,
    smois,
    ust,
    isltyp,
    vegfra,
    dt,
    emis_dust,
    soil_top_cat,
    erod,
    sf_surface_physics=2,
    imod=3,
    psd=None,
    settings=None,
):
    """Port U:22-271 on WRF arrays; return updated state and diagnostics.

    `dust` is [x,level,y,5] in ug/kg; `rho_phy` and `dz8w` are [x,level,y];
    `smois` and `soil_top_cat` are [x,soil,y]. Surface arrays are [x,y],
    `erod` is [x,y,3], and `emis_dust` is [x,y,5] in kg m-2.
    """
    settings = Shao2011Settings() if settings is None else settings
    shape = np.shape(ust)
    if (
        len(shape) != 2
        or np.ndim(dust) != 4
        or np.shape(dust)[1] < 1
        or np.shape(dust) != (*shape[:1], np.shape(dust)[1], *shape[1:], N_BINS)
        or np.shape(emis_dust) != (*shape, N_BINS)
        or np.ndim(rho_phy) != 3
        or np.ndim(dz8w) != 3
        or np.shape(rho_phy)[1] < 1
        or np.shape(dz8w)[1] < 1
        or np.shape(rho_phy)[0::2] != shape
        or np.shape(dz8w)[0::2] != shape
        or np.ndim(smois) != 3
        or np.shape(smois)[1] < 1
        or np.shape(smois)[0::2] != shape
        or np.shape(soil_top_cat)
        not in ((*shape[:1], 12, *shape[1:]), (*shape[:1], 16, *shape[1:]))
    ):
        raise ValueError("[dust_shao2011] inconsistent WRF driver shapes")
    if psd is None:
        psd = initialize(settings=settings)
    soil = np.moveaxis(np.asarray(soil_top_cat)[:, :12, :], 1, -1)
    wet = np.array(smois, dtype=float, copy=True)
    flux, capped, ust_min, rough, moist = _flux_and_diagnostics(
        ust,
        wet[:, 0, :],
        vegfra,
        erod,
        soil,
        np.asarray(rho_phy)[:, 0, :],
        psd,
        sf_surface_physics,
        imod,
        isltyp,
        settings,
    )
    active = np.asarray(erod).sum(axis=-1) > 0.0
    wet[:, 0, :] = capped
    updated = np.array(dust, dtype=float, copy=True)
    emitted = np.array(emis_dust, dtype=float, copy=True)
    mass = flux * dt
    updated[:, 0, :, :][active] = (
        updated[:, 0, :, :][active] * 1.0e-9
        + mass[active]
        / np.asarray(dz8w)[:, 0, :][active, None]
        / np.asarray(rho_phy)[:, 0, :][active, None]
    ) * 1.0e9
    emitted[active] += mass[active]
    emitted[~active] = 0.0
    return dict(
        dust=updated,
        emis_dust=emitted,
        smois=wet,
        ust_t_min=ust_min,
        rough_cor=rough,
        smois_cor=moist,
    )


# ===========================================================================
# Box-model demo
# Run with: python dust_shao2011.py
# ===========================================================================
