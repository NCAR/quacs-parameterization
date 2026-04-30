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
dependencies beyond NumPy. Runs as a script (python dust_ginoux.py)
or as a pytest test suite (pytest dust_ginoux.py -v).

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
import pytest

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80616    # m/s^2
LAND = 1.0        # oro value identifying land cells

# ---------------------------------------------------------------------------
# Bin parameters — source: DU2G_instance.F90
# ---------------------------------------------------------------------------

# Effective particle radius per bin [m]
# Source: DU2G_instance.F90, L26  particle_radius_microns = (0.73,1.4,2.4,4.5,8.0)
RADIUS_DEFAULT = np.array([0.73, 1.4, 2.4, 4.5, 8.0]) * 1.0e-6   # [m]

# Bin radius boundaries [m]
# Source: DU2G_instance.F90, L28-30
RADIUS_LOWER = np.array([0.1,  1.0, 1.8, 3.0,  6.0]) * 1.0e-6    # [m]
RADIUS_UPPER = np.array([1.0,  1.8, 3.0, 6.0, 10.0]) * 1.0e-6    # [m]

# Particle bulk density per bin [kg/m^3]
# Source: DU2G_instance.F90, L33  particle_density = (2500,2650,2650,2650,2650)
RHOP_DEFAULT = np.array([2500., 2650., 2650., 2650., 2650.])       # [kg/m^3]

# Per-bin scaling factor [-]
# Source: DU2G_instance.F90, L43  source_fraction = (0.1,0.25,0.25,0.25,0.25)
# Note: sums to 1.1, not 1.0. In GOCART2G, sfra is not used inside the
# emission formula itself — per-bin variation comes from du_src + ipoint.
# SP_DEFAULT reproduces the Fortran values exactly as a reference baseline.
# See the module docstring for a full explanation.
SP_DEFAULT = np.array([0.1, 0.25, 0.25, 0.25, 0.25])              # sum = 1.1

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
    d  = (2.0 * radius)[np.newaxis, np.newaxis, :]    # diameter [m], [1, 1, n_bins]
    rp = rhop[np.newaxis, np.newaxis, :]               # [1, 1, n_bins]
    ra = rhoa[:, :, np.newaxis]                        # [x, y, 1]

    u_t = (
        0.13
        * np.sqrt(rp * GRAV * d / ra)
        * np.sqrt(1.0 + 6.0e-7 / (rp * GRAV * d**2.5))
        / np.sqrt(1.928 * (1331.0 * (100.0 * d)**1.56 + 0.38)**0.092 - 1.0)
    )
    return u_t    # [x, y, n_bins]

# ---------------------------------------------------------------------------
# Initialization: load static external datasets once at model start
# ---------------------------------------------------------------------------
def initialize(erodibility_file, sp_file, oro_file, frlake_file,
               n_bins=len(SP_DEFAULT)):
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
    files  = [erodibility_file, sp_file, oro_file, frlake_file]
    labels = ['erodibility_file', 'sp_file', 'oro_file', 'frlake_file']
    for path, label in zip(files, labels):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[dust_ginoux] {label} not found: {path}")

    S      = np.load(erodibility_file)
    s_p    = np.load(sp_file)
    oro    = np.load(oro_file)
    frlake = np.load(frlake_file)

    if S.ndim != 2:
        raise ValueError(f"[dust_ginoux] S must be 2D [x,y], got {S.shape}")
    if not (S.min() >= 0.0 and S.max() <= 1.0):
        raise ValueError(f"[dust_ginoux] S must be in [0,1]; got [{S.min():.3f}, {S.max():.3f}]")
    if s_p.shape != (n_bins,):
        raise ValueError(f"[dust_ginoux] s_p must have shape ({n_bins},), got {s_p.shape}")
    if (s_p < 0.0).any():
        raise ValueError("[dust_ginoux] s_p contains negative values")
    if oro.shape != S.shape:
        raise ValueError(f"[dust_ginoux] oro shape {oro.shape} != S shape {S.shape}")
    if frlake.shape != S.shape:
        raise ValueError(f"[dust_ginoux] frlake shape {frlake.shape} != S shape {S.shape}")
    if not (frlake.min() >= 0.0 and frlake.max() <= 1.0):
        raise ValueError("[dust_ginoux] frlake must be in [0, 1]")

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
    u_thresh     = np.maximum(u_t * (1.2 + 0.2 * np.log10(gwet_clamped)), 0.0)

    # --- combined emission mask: land, not too wet, wind above threshold -----
    emit_mask   = (oro == LAND) & (gwet < 0.5) & (w10m > u_thresh)

    # --- wind excess and lake fraction correction ----------------------------
    wind_excess = np.where(emit_mask, w10m - u_thresh, 0.0)   # [m/s]
    lake_factor = 1.0 - frlake                                 # [-]

    # --- emission flux: [x, y, 1] * [1, 1, n_bins] -> [x, y, n_bins] -------
    base_flux = C * S * lake_factor * w10m**2 * wind_excess    # [kg m^-2 s^-1]
    F_p       = base_flux[:, :, np.newaxis] * s_p[np.newaxis, np.newaxis, :]

    return F_p    # [kg m^-2 s^-1]

# ---------------------------------------------------------------------------
# Pytest test suite
# Run with: pytest dust_ginoux.py -v
# ---------------------------------------------------------------------------

@pytest.fixture
def cell():
    """Single land cell, moderate moisture, wind above threshold."""
    return dict(
        w10m   = np.array([[12.0]]),
        u_t    = np.array([[ 6.0]]),
        gwet   = np.array([[ 0.1]]),
        oro    = np.array([[ LAND]]),
        frlake = np.array([[ 0.0]]),
        S      = np.array([[ 0.5]]),
        s_p    = SP_DEFAULT.copy(),
        C      = C_DEFAULT,
    )

# --- shape and basic contract ---

def test_output_shape(cell):
    """F_p must have shape [x, y, n_bins]."""
    assert dust_emission(**cell).shape == (1, 1, 5)

def test_all_fluxes_nonnegative(cell):
    """No negative emission fluxes anywhere."""
    assert (dust_emission(**cell) >= 0.0).all()

def test_bin_fluxes_scale_with_sp(cell):
    """
    Each bin's flux must be proportional to its s_p value.
    F_p[bin_i] / F_p[bin_j] == s_p[i] / s_p[j] for all pairs.
    This holds regardless of whether s_p sums to 1.
    """
    F_p = dust_emission(**cell)[0, 0, :]
    for i in range(len(F_p) - 1):
        ratio_flux = F_p[i] / F_p[i+1]
        ratio_sp   = cell['s_p'][i] / cell['s_p'][i+1]
        assert np.isclose(ratio_flux, ratio_sp, rtol=1e-10), \
            f"Bin {i}/{i+1} flux ratio {ratio_flux:.4f} != s_p ratio {ratio_sp:.4f}"

def test_sp_default_matches_fortran():
    """
    SP_DEFAULT must match DU2G_instance.F90, L43 exactly.
    Note: source_fraction sums to 1.1 in the Fortran; this is expected.
    See module docstring for explanation of why s_p is not a strict partition.
    """
    expected = np.array([0.1, 0.25, 0.25, 0.25, 0.25])
    assert np.allclose(SP_DEFAULT, expected), \
        f"SP_DEFAULT {SP_DEFAULT} does not match Fortran source_fraction {expected}"
    assert np.isclose(SP_DEFAULT.sum(), 1.1), \
        f"Expected sum 1.1 (as in Fortran); got {SP_DEFAULT.sum():.4f}"

# --- analytic correctness ---

def test_total_flux_matches_analytic(cell):
    """Total flux must match the analytic formula exactly."""
    c            = cell
    gwet_clamped = max(c['gwet'][0,0], 1e-3)
    u_thresh     = max(c['u_t'][0,0] * (1.2 + 0.2 * np.log10(gwet_clamped)), 0.0)
    # base_flux * sum(s_p) = total across bins
    base_flux    = c['C'] * c['S'][0,0] * c['w10m'][0,0]**2 * (c['w10m'][0,0] - u_thresh)
    expected     = base_flux * c['s_p'].sum()
    F_p          = dust_emission(**c)
    assert np.isclose(F_p[0, 0, :].sum(), expected, rtol=1e-10)

# --- emission cutoffs: land mask ---

def test_zero_flux_over_ocean():
    """Emission must be zero over non-land cells."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]), u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),  oro=np.array([[2.0]]),
        frlake=np.array([[0.0]]), S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()

# --- emission cutoffs: soil moisture ---

def test_zero_flux_wet_soil():
    """Emission must be zero when gwet >= 0.5."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]), u_t=np.array([[5.0]]),
        gwet=np.array([[0.5]]),  oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]), S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()

def test_moisture_correction_raises_threshold():
    """Moister soil raises effective threshold, reducing flux."""
    base = dict(w10m=np.array([[10.0]]), u_t=np.array([[6.0]]),
                oro=np.array([[LAND]]),  frlake=np.array([[0.0]]),
                S=np.array([[0.5]]),     s_p=SP_DEFAULT)
    F_dry   = dust_emission(gwet=np.array([[0.01]]), **base)
    F_moist = dust_emission(gwet=np.array([[0.40]]), **base)
    assert F_moist[0,0,:].sum() < F_dry[0,0,:].sum()

# --- emission cutoffs: wind threshold ---

def test_zero_flux_below_threshold():
    """Wind below adjusted threshold must give zero flux."""
    F_p = dust_emission(
        w10m=np.array([[3.0]]), u_t=np.array([[8.0]]),
        gwet=np.array([[0.1]]), oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]), S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()

# --- lake fraction ---

def test_full_lake_gives_zero_flux():
    """frlake = 1.0 must give zero flux."""
    F_p = dust_emission(
        w10m=np.array([[12.0]]), u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),  oro=np.array([[LAND]]),
        frlake=np.array([[1.0]]), S=np.array([[0.5]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()

def test_lake_fraction_scales_flux():
    """frlake = 0.5 must give exactly half the flux of frlake = 0."""
    base = dict(w10m=np.array([[12.0]]), u_t=np.array([[5.0]]),
                gwet=np.array([[0.1]]),  oro=np.array([[LAND]]),
                S=np.array([[0.5]]),     s_p=SP_DEFAULT)
    F_no_lake   = dust_emission(frlake=np.array([[0.0]]), **base)
    F_half_lake = dust_emission(frlake=np.array([[0.5]]), **base)
    assert np.isclose(F_half_lake.sum(), F_no_lake.sum() * 0.5, rtol=1e-10)

# --- zero erodibility ---

def test_zero_erodibility_gives_zero_flux():
    """S = 0 must give zero flux regardless of all other conditions."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]), u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),  oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]), S=np.array([[0.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()

# --- scaling laws ---

def test_flux_scales_linearly_with_C(cell):
    """Doubling C must exactly double the total flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, 'C': cell['C'] * 2.0})
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-10)

def test_flux_scales_linearly_with_S(cell):
    """Doubling S must exactly double the total flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, 'S': cell['S'] * 2.0})
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-10)

# --- vectorized grid ---

def test_vectorized_grid():
    """Must run correctly on a full 2D grid with mixed land/ocean cells."""
    nx, ny, nb = 12, 10, 5
    rng = np.random.default_rng(0)
    F_p = dust_emission(
        w10m   = rng.uniform(0,    20, (nx, ny)),
        u_t    = rng.uniform(4,     8, (nx, ny)),
        gwet   = rng.uniform(0,  0.45, (nx, ny)),
        oro    = rng.choice([LAND, 2.0], (nx, ny)),
        frlake = rng.uniform(0,     1, (nx, ny)),
        S      = rng.uniform(0,     1, (nx, ny)),
        s_p    = SP_DEFAULT,
    )
    assert F_p.shape == (nx, ny, nb)
    assert (F_p >= 0.0).all()

# --- threshold_velocity ---

def test_threshold_velocity_shape():
    """threshold_velocity must return shape [x, y, n_bins]."""
    u_t = threshold_velocity(RADIUS_DEFAULT, RHOP_DEFAULT, np.ones((4, 6)) * 1.2)
    assert u_t.shape == (4, 6, 5)

def test_threshold_velocity_positive():
    """All threshold velocities must be positive for physical inputs."""
    u_t = threshold_velocity(RADIUS_DEFAULT, RHOP_DEFAULT, np.ones((3, 3)) * 1.2)
    assert (u_t > 0.0).all()

def test_threshold_decreases_with_particle_size():
    """
    In the fine-dust regime (0.73-8 um), threshold decreases with particle
    size. Cohesive forces dominate small particles (Marticorena 1995).
    The U-curve minimum is near ~100 um; all GOCART bins are well below it.
    """
    u_t = threshold_velocity(RADIUS_DEFAULT, RHOP_DEFAULT, np.ones((1, 1)) * 1.2)
    assert np.all(np.diff(u_t[0, 0, :]) < 0), \
        f"Expected decreasing threshold across bins; got {u_t[0,0,:]}"

# --- source-verified parameter guards ---

def test_radius_bin1_is_073_micron():
    """
    Bin 1 radius must be 0.73 um (DU2G_instance.F90, L26).
    Guards against reverting to the earlier placeholder value of 0.5 um.
    """
    assert np.isclose(RADIUS_DEFAULT[0], 0.73e-6), \
        f"Bin 1 radius should be 0.73e-6 m; got {RADIUS_DEFAULT[0]:.2e}"

def test_particle_density_per_bin():
    """
    Bin 1 density must be 2500 kg/m^3; bins 2-5 must be 2650 kg/m^3.
    Source: DU2G_instance.F90, L33.
    """
    assert RHOP_DEFAULT[0] == 2500.0, f"Bin 1 density: expected 2500, got {RHOP_DEFAULT[0]}"
    assert np.all(RHOP_DEFAULT[1:] == 2650.0), \
        f"Bins 2-5 density: expected 2650, got {RHOP_DEFAULT[1:]}"

# --- initialize() ---

def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("no_S.npy", "no_sp.npy", "no_oro.npy", "no_frlake.npy")

def test_initialize_sp_wrong_shape(tmp_path):
    """initialize() must raise ValueError if s_p has wrong shape."""
    S_f  = tmp_path / "S.npy"
    sp_f = tmp_path / "sp.npy"
    or_f = tmp_path / "oro.npy"
    fl_f = tmp_path / "frl.npy"
    np.save(S_f,  np.ones((4,4)) * 0.5)
    np.save(sp_f, np.array([0.25, 0.25, 0.25, 0.25]))   # 4 bins, expected 5
    np.save(or_f, np.ones((4,4)))
    np.save(fl_f, np.zeros((4,4)))
    with pytest.raises(ValueError, match="shape"):
        initialize(str(S_f), str(sp_f), str(or_f), str(fl_f))

def test_initialize_valid(tmp_path):
    """initialize() must return correct shapes for valid inputs."""
    S_f  = tmp_path / "S.npy"
    sp_f = tmp_path / "sp.npy"
    or_f = tmp_path / "oro.npy"
    fl_f = tmp_path / "frl.npy"
    np.save(S_f,  np.ones((6,6)) * 0.4)
    np.save(sp_f, SP_DEFAULT)
    np.save(or_f, np.ones((6,6)))
    np.save(fl_f, np.zeros((6,6)))
    S, s_p, oro, frlake = initialize(str(S_f), str(sp_f), str(or_f), str(fl_f))
    assert S.shape      == (6, 6)
    assert s_p.shape    == (5,)
    assert oro.shape    == (6, 6)
    assert frlake.shape == (6, 6)

# ---------------------------------------------------------------------------
# Box-model demo
# Run with: python dust_ginoux.py
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # synthetic single-cell inputs
    w10m_val = 12.0    # m/s   10-m wind speed
    gwet_val =  0.1    # [-]   gravimetric soil moisture
    S_val    =  0.5    # [-]   erodibility
    rhoa_val =  1.2    # kg/m^3 near-surface air density

    # compute dry-soil threshold via Marticorena & Bergametti (1995)
    u_t_bins = threshold_velocity(
        RADIUS_DEFAULT, RHOP_DEFAULT, np.array([[rhoa_val]])
    )

    print("Bin parameters (source: DU2G_instance.F90):")
    print(f"  Radii [um]:        {RADIUS_DEFAULT * 1e6}")
    print(f"  Densities [kg/m3]: {RHOP_DEFAULT}")
    print(f"  s_p [-]:           {SP_DEFAULT}  (sum={SP_DEFAULT.sum():.2f})")
    print(f"  u_t [m/s]:         {u_t_bins[0,0,:].round(3)}")
    print()

    # use bin-mean threshold for single-cell demo
    u_t_mean = float(u_t_bins[0, 0, :].mean())

    F_p = dust_emission(
        w10m   = np.array([[w10m_val]]),
        u_t    = np.array([[u_t_mean]]),
        gwet   = np.array([[gwet_val]]),
        oro    = np.array([[LAND]]),
        frlake = np.array([[0.0]]),
        S      = np.array([[S_val]]),
        s_p    = SP_DEFAULT,
    )

    # flux-to-tendency time loop (QUACS doc v1, Section 7A)
    rho_air = 1.2    # kg/m^3
    dz      = 50.0   # m
    dt      = 300.0  # s
    x       = np.zeros(len(SP_DEFAULT))

    print("QUACS dust emission (GOCART) — box model demo")
    print("=" * 56)
    print(f"  w10m={w10m_val} m/s | gwet={gwet_val} | S={S_val} | C={C_DEFAULT}")
    print(f"  u_t (bin-mean) = {u_t_mean:.3f} m/s")
    print()
    print(f"{'Step':>5}  {'F_total [kg/m2/s]':>20}  {'x_total [kg/kg]':>18}")
    print("-" * 50)
    for step in range(1, 11):
        tendency = F_p[0, 0, :] / (rho_air * dz)
        x       += tendency * dt
        print(f"{step:>5}  {F_p[0,0,:].sum():>20.4e}  {x.sum():>18.4e}")

    print()
    print("Run tests with: pytest dust_ginoux.py -v")