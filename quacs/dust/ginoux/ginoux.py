"""
Dust Emission Scheme: Ginoux et al. (2001) / GOCART
====================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Ginoux, P., et al. (2001). Sources and distributions of dust aerosols
    simulated with the GOCART model. J. Geophys. Res., 106, 20255-20273.
    https://doi.org/10.1029/2000JD000053

    Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric
    dust cycle: 1. Design of a soil-derived dust emission scheme.
    J. Geophys. Res., 100, 16415-16430.
    https://doi.org/10.1029/95JD00690

    Fortran reference: MPAS-GOCART2G
      DU2G_GridCompMod.F90  subroutine DustEmissionGOCART2G_revised (L728-L825)
      DU2G_instance.F90     all bin parameters (radii, densities, source_fraction)
    https://github.com/PACE-DAAQ/MPAS-GOCART2G

QUACS Contract
--------------
Standalone Python module. No MPAS-A, no config files, no external
dependencies beyond NumPy. The tests are in tests/test_ginoux.py. The
MUSICA box model is in box_model.py.

Pre-API note:
    Once the QUACS scripting API is finalized, the expected signature is:
        @scheme(category='dust_emission')
        def dust_emission_gocart(state: ModelState, config: Config) -> SurfaceFlux:
    Until then, plain NumPy arrays are used. This script is the scientific
    specification that drives that API design.

Physics summary (Ginoux et al. 2001, Eq. 1):
    F_p = C * S * s_p * (1 - frlake) * w10m^2 * max(w10m - u_thresh, 0)

    where u_thresh is the moisture-adjusted threshold (Ginoux 2001):
        u_thresh = max(0, u_t * (1.2 + 0.2 * log10(max(gwet, 1e-3))))

    Emission is zeroed over ocean/sea-ice (oro != 1), over lakes (frlake),
    and when soil moisture is too high (gwet >= 0.5).

    The dry-soil threshold u_t is computed per bin from particle physics
    via threshold_velocity() following Marticorena & Bergametti (1995).

Inputs (dynamic, per time step from MPAS-A):
    w10m  : np.ndarray [x, y]     10-m wind speed magnitude [m/s]
    gwet  : np.ndarray [x, y]     gravimetric soil moisture [-], in [0, 1]

Inputs (static, loaded once at initialization):
    S      : np.ndarray [x, y]     source erodibility factor [-], in [0, 1]
    s_p    : np.ndarray [n_bins]   scaling factor per size bin [-]
    oro    : np.ndarray [x, y]     land-water mask (1.0=land, 2.0=ocean, etc.)
    frlake : np.ndarray [x, y]     lake fraction [-], in [0, 1]
    radius : np.ndarray [n_bins]   effective particle radius per bin [m]
    rhop   : np.ndarray [n_bins]   particle density per bin [kg/m^3]
    rhoa   : np.ndarray [x, y]     near-surface air density [kg/m^3]

Configurable:
    C  : float   empirical scaling constant.
                 In GOCART2G (DU2G_instance.F90, L42) Ch_DU is resolution-
                 dependent: [0.3, 0.3, 0.11, 0.11, 0.11, 0.088] for six
                 grid resolutions (a through f). C_DEFAULT uses the
                 coarsest-resolution value (0.088) as a standalone default.
                 empirical scaling constant [kg s^2 m^-5]
                                default: 1.0e-9 (Ginoux 2001)
                 Select the appropriate value for your target grid resolution.

Output:
    F_p : np.ndarray [x, y, n_bins]   surface emission flux [kg m^-2 s^-1]
                                      QUACS runtime converts to tendency:
                                      dx/dt = F_p / (rho_air * dz)

Note on s_p vs. source_fraction:
    In GOCART2G, per-bin variation in emissions comes from du_src [x, y, n_src]
    combined with ipoint [n_bins] (DU2G_instance.F90, L46: ipoint = [3,2,2,2,2]),
    which maps each bin to an erodibility source category. The variable sfra
    (source_fraction) is passed into DustEmissionGOCART2G_revised but is NOT
    used inside the emission formula itself.

    This script simplifies that mechanism into a single s_p [n_bins] multiplier.
    SP_DEFAULT reproduces the Fortran source_fraction values directly; note that
    these sum to 1.1, not 1.0, consistent with the Fortran source. This is not
    a normalization error — s_p here is a per-bin scaling factor, not strictly
    a fractional partition.
"""

import os

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80616  # m/s^2
LAND = 1.0  # oro value identifying land cells

# ---------------------------------------------------------------------------
# Bin parameters — source: DU2G_instance.F90
# ---------------------------------------------------------------------------

# Effective particle radius per bin [m]
# Source: DU2G_instance.F90, L26  particle_radius_microns = (0.73,1.4,2.4,4.5,8.0)
RADIUS_DEFAULT = np.array([0.73, 1.4, 2.4, 4.5, 8.0]) * 1.0e-6  # [m]

# Bin radius boundaries [m]
# Source: DU2G_instance.F90, L28-30
RADIUS_LOWER = np.array([0.1, 1.0, 1.8, 3.0, 6.0]) * 1.0e-6  # [m]
RADIUS_UPPER = np.array([1.0, 1.8, 3.0, 6.0, 10.0]) * 1.0e-6  # [m]

# Particle bulk density per bin [kg/m^3]
# Source: DU2G_instance.F90, L33  particle_density = (2500,2650,2650,2650,2650)
RHOP_DEFAULT = np.array([2500.0, 2650.0, 2650.0, 2650.0, 2650.0])  # [kg/m^3]

# Per-bin scaling factor [-]
# Source: DU2G_instance.F90, L43  source_fraction = (0.1,0.25,0.25,0.25,0.25)
# Note: sums to 1.1, not 1.0. In GOCART2G, sfra is not used inside the
# emission formula itself — per-bin variation comes from du_src + ipoint.
# SP_DEFAULT reproduces the Fortran values exactly as a reference baseline.
# See the module docstring for a full explanation.
SP_DEFAULT = np.array([0.1, 0.25, 0.25, 0.25, 0.25])  # sum = 1.1

# Empirical scaling constant
# Source: DU2G_instance.F90, L42
#   Ch_DU = (0.3, 0.3, 0.11, 0.11, 0.11, 0.088) for resolutions a-f.
# C_DEFAULT uses the coarsest-resolution value as a conservative standalone default.
# Previous original Ginoux (2001) value was 1.0e-9 with different units.
C_DEFAULT = 0.088


# ---------------------------------------------------------------------------
# Threshold velocity: Marticorena & Bergametti (1995)
# ---------------------------------------------------------------------------
def threshold_velocity(radius, rhop, rhoa):
    """
    Compute the dry-soil threshold wind speed per size bin and grid cell.

    Follows Marticorena & Bergametti (1995), as implemented in GOCART2G
    DustEmissionGOCART2G_revised (DU2G_GridCompMod.F90, L800-L802).

    In the fine-dust range (0.73-8 um), threshold DECREASES with particle
    size. Cohesive inter-particle forces dominate small particles, making
    them harder to lift. The classic Marticorena U-curve minimum is near
    ~100 um (sand); all GOCART bins are well below that inflection point.

    Parameters
    ----------
    radius : np.ndarray [n_bins]   effective particle radius per bin [m]
    rhop   : np.ndarray [n_bins]   particle density per bin [kg/m^3]
    rhoa   : np.ndarray [x, y]     near-surface air density [kg/m^3]

    Returns
    -------
    u_t : np.ndarray [x, y, n_bins]   dry-soil threshold wind speed [m/s]
    """
    d = (2.0 * radius)[np.newaxis, np.newaxis, :]  # diameter [m], [1, 1, n_bins]
    rp = rhop[np.newaxis, np.newaxis, :]  # [1, 1, n_bins]
    ra = rhoa[:, :, np.newaxis]  # [x, y, 1]

    u_t = (
        0.13
        * np.sqrt(rp * GRAV * d / ra)
        * np.sqrt(1.0 + 6.0e-7 / (rp * GRAV * d**2.5))
        / np.sqrt(1.928 * (1331.0 * (100.0 * d) ** 1.56 + 0.38) ** 0.092 - 1.0)
    )
    return u_t  # [x, y, n_bins]


# ---------------------------------------------------------------------------
# Initialization: load static external datasets once at model start
# ---------------------------------------------------------------------------
def initialize(erodibility_file, sp_file, oro_file, frlake_file, n_bins=len(SP_DEFAULT)):
    """
    Load and validate static input fields. Called once before the time loop.

    Parameters
    ----------
    erodibility_file : str   Path to S [x, y], values in [0, 1].
    sp_file          : str   Path to s_p [n_bins], per-bin scaling factors.
    oro_file         : str   Path to land-water mask [x, y] (1.0 = land).
    frlake_file      : str   Path to lake fraction [x, y], values in [0, 1].
    n_bins           : int   Expected number of size bins.

    Returns
    -------
    S      : np.ndarray [x, y]
    s_p    : np.ndarray [n_bins]
    oro    : np.ndarray [x, y]
    frlake : np.ndarray [x, y]
    """
    files = [erodibility_file, sp_file, oro_file, frlake_file]
    labels = ["erodibility_file", "sp_file", "oro_file", "frlake_file"]
    for path, label in zip(files, labels):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[ginoux] {label} not found: {path}")

    S = np.load(erodibility_file)
    s_p = np.load(sp_file)
    oro = np.load(oro_file)
    frlake = np.load(frlake_file)

    if S.ndim != 2:
        raise ValueError(f"[ginoux] S must be 2D [x,y], got {S.shape}")
    if not (S.min() >= 0.0 and S.max() <= 1.0):
        raise ValueError(f"[ginoux] S must be in [0,1]; got [{S.min():.3f}, {S.max():.3f}]")
    if s_p.shape != (n_bins,):
        raise ValueError(f"[ginoux] s_p must have shape ({n_bins},), got {s_p.shape}")
    if (s_p < 0.0).any():
        raise ValueError("[ginoux] s_p contains negative values")
    if oro.shape != S.shape:
        raise ValueError(f"[ginoux] oro shape {oro.shape} != S shape {S.shape}")
    if frlake.shape != S.shape:
        raise ValueError(f"[ginoux] frlake shape {frlake.shape} != S shape {S.shape}")
    if not (frlake.min() >= 0.0 and frlake.max() <= 1.0):
        raise ValueError("[ginoux] frlake must be in [0, 1]")

    return S, s_p, oro, frlake


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def dust_emission(w10m, u_t, gwet, oro, frlake, S, s_p, C=C_DEFAULT):
    """
    Compute dust emission flux: Ginoux et al. (2001), Eq. 1.

    Physics (GOCART2G / DustEmissionGOCART2G_revised):

        u_thresh = max(0, u_t * (1.2 + 0.2 * log10(max(gwet, 1e-3))))
        F_p = C * S * (1 - frlake) * s_p * w10m^2 * max(w10m - u_thresh, 0)

    Emission is zero over non-land cells (oro != 1.0) and when gwet >= 0.5.

    Parameters
    ----------
    w10m   : np.ndarray [x, y]            10-m wind speed [m/s]
    u_t    : np.ndarray [x, y] or scalar  dry-soil threshold wind speed [m/s],
                                          typically from threshold_velocity()
    gwet   : np.ndarray [x, y]            gravimetric soil moisture [-]
    oro    : np.ndarray [x, y]            land-water mask (1.0 = land)
    frlake : np.ndarray [x, y]            lake fraction [-]
    S      : np.ndarray [x, y]            source erodibility factor [-]
    s_p    : np.ndarray [n_bins]          per-bin scaling factor [-]
    C      : float                        scaling constant (see C_DEFAULT note)

    Returns
    -------
    F_p : np.ndarray [x, y, n_bins]   dust emission flux [kg m^-2 s^-1]
    """
    # --- soil moisture adjusted threshold (Ginoux 2001) ----------------------
    gwet_clamped = np.maximum(gwet, 1.0e-3)
    u_thresh = np.maximum(u_t * (1.2 + 0.2 * np.log10(gwet_clamped)), 0.0)

    # --- combined emission mask: land, not too wet, wind above threshold -----
    emit_mask = (oro == LAND) & (gwet < 0.5) & (w10m > u_thresh)

    # --- wind excess and lake fraction correction ----------------------------
    wind_excess = np.where(emit_mask, w10m - u_thresh, 0.0)  # [m/s]
    lake_factor = 1.0 - frlake  # [-]

    # --- emission flux: [x, y, 1] * [1, 1, n_bins] -> [x, y, n_bins] -------
    base_flux = C * S * lake_factor * w10m**2 * wind_excess  # [kg m^-2 s^-1]
    F_p = base_flux[:, :, np.newaxis] * s_p[np.newaxis, np.newaxis, :]

    return F_p  # [kg m^-2 s^-1]
