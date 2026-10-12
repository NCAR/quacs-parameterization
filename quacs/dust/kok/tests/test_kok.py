"""Unit tests for the kok dust emission scheme.

These tests came from the original standalone module by Sue Park.
"""

import numpy as np
import pytest

from quacs.dust.kok.kok import (
    C_KAPPA,
    C_TUNE_DEFAULT,
    CD0,
    CE,
    DP_DEFAULT,
    FECAN_A_DEFAULT,
    GRAV,
    KAPPA_MAX,
    LAND,
    RHO_A0,
    RHOP_DEFAULT,
    SL_A,
    SL_GAMMA,
    U_ST0,
    VAI_THR,
    bare_fraction,
    dust_emission,
    erodibility_coefficient,
    fragmentation_exponent,
    initialize,
    moisture_factor,
    standardized_threshold,
    threshold_velocity_dry,
)


@pytest.fixture
def cell():
    """Single land cell: moderate moisture, wind above threshold."""
    return dict(
        ustar=np.array([[0.5]]),  # m/s  — well above typical u_ft0
        w=np.array([[0.02]]),  # kg/kg — dry, f_m ~ 1
        rho_a=np.array([[1.2]]),  # kg/m^3
        vai=np.array([[0.05]]),  # [-]  — sparse veg, f_bare > 0
        snowdp=np.array([[0.0]]),  # m    — no snow
        f_clay=np.array([[0.15]]),  # [-]  — typical desert soil
        oro=np.array([[LAND]]),
        Dp=DP_DEFAULT,
        rhop=RHOP_DEFAULT,
        a=FECAN_A_DEFAULT,
        vai_thr=VAI_THR,
        C_tune=C_TUNE_DEFAULT,
    )


# --- shape and basic contract ---


def test_output_shape(cell):
    """F_d must have shape [x, y]."""
    F_d = dust_emission(**cell)
    assert F_d.shape == (1, 1)


def test_flux_nonnegative(cell):
    """Emission flux must be non-negative everywhere."""
    assert (dust_emission(**cell) >= 0.0).all()


def test_flux_positive_for_emitting_cell(cell):
    """A land cell with high wind, dry soil, sparse veg must emit dust."""
    assert dust_emission(**cell)[0, 0] > 0.0


# --- analytic building-block checks ---


def test_threshold_velocity_dry_positive():
    """Dry threshold must be positive for all physical inputs."""
    rho_a = np.ones((3, 4)) * 1.2
    u_ft0 = threshold_velocity_dry(Dp=DP_DEFAULT, rhop=RHOP_DEFAULT, rho_a=rho_a)
    assert (u_ft0 > 0.0).all()


def test_threshold_velocity_dry_scalar():
    """threshold_velocity_dry must accept scalar rho_a (None → RHO_A0)."""
    u_ft0 = threshold_velocity_dry()
    assert float(u_ft0) > 0.0


def test_moisture_factor_dry_soil():
    """For very dry soil (w << w_t), f_m must equal 1.0."""
    w = np.array([[0.001]])
    f_clay = np.array([[0.15]])
    f_m = moisture_factor(w, f_clay)
    assert np.isclose(f_m[0, 0], 1.0, atol=1e-10)


def test_moisture_factor_wet_soil_gt_1():
    """For moist soil (w > w_t), f_m must be > 1."""
    f_clay = np.array([[0.15]])
    # w_t ≈ 0.01 * 1 * (17*0.15 + 14*0.15^2) ≈ 0.0287
    w = np.array([[0.10]])  # well above w_t
    f_m = moisture_factor(w, f_clay)
    assert f_m[0, 0] > 1.0


def test_moisture_factor_increases_with_w():
    """f_m must increase monotonically with soil moisture above w_t."""
    f_clay = np.array([[0.15]])
    w_vals = np.array([0.05, 0.10, 0.20, 0.40])
    f_m_vals = np.array([moisture_factor(np.array([[w]]), f_clay)[0, 0] for w in w_vals])
    assert np.all(np.diff(f_m_vals) >= 0.0)


def test_standardized_threshold_scales_with_rho():
    """u_st must increase when rho_a increases (denser air → higher u_st)."""
    u_ft = np.array([[0.3]])
    rho_1 = np.array([[1.0]])
    rho_2 = np.array([[1.4]])
    assert standardized_threshold(u_ft, rho_2) > standardized_threshold(u_ft, rho_1)


def test_erodibility_decreases_with_u_st():
    """C_d must decrease as u_st increases (less erodible at higher threshold)."""
    u_st_low = np.array([[0.16]])
    u_st_high = np.array([[0.40]])
    assert erodibility_coefficient(u_st_low) > erodibility_coefficient(u_st_high)


def test_erodibility_at_u_st0():
    """At u_st = u_st0, C_d must equal C_d0 exactly."""
    u_st = np.array([[U_ST0]])
    assert np.isclose(erodibility_coefficient(u_st)[0, 0], CD0, rtol=1e-12)


def test_fragmentation_exponent_capped():
    """kappa must never exceed KAPPA_MAX."""
    u_st = np.linspace(U_ST0, 1.5, 100).reshape(10, 10)
    kappa = fragmentation_exponent(u_st)
    assert (kappa <= KAPPA_MAX + 1e-12).all()


def test_fragmentation_exponent_at_u_st0():
    """At u_st = u_st0, kappa must equal 0."""
    u_st = np.array([[U_ST0]])
    assert np.isclose(fragmentation_exponent(u_st)[0, 0], 0.0, atol=1e-12)


def test_bare_fraction_no_veg():
    """With VAI = 0, f_bare must be 1."""
    assert np.isclose(bare_fraction(np.array([[0.0]]))[0, 0], 1.0)


def test_bare_fraction_at_threshold():
    """With VAI = VAI_THR, f_bare must be 0."""
    assert np.isclose(bare_fraction(np.array([[VAI_THR]]))[0, 0], 0.0, atol=1e-12)


def test_bare_fraction_above_threshold():
    """With VAI > VAI_THR, f_bare must be 0 (clipped)."""
    assert bare_fraction(np.array([[VAI_THR * 2]]))[0, 0] == 0.0


# --- emission cutoffs ---


def test_zero_flux_over_ocean(cell):
    """Emission must be zero over non-land cells."""
    c = {**cell, "oro": np.array([[2.0]])}
    assert dust_emission(**c)[0, 0] == 0.0


def test_zero_flux_snow_covered(cell):
    """Emission must be zero when snowdp > SNOW_THRESH."""
    c = {**cell, "snowdp": np.array([[0.05]])}
    assert dust_emission(**c)[0, 0] == 0.0


def test_zero_flux_below_wind_threshold(cell):
    """Emission must be zero when ustar <= u_ft."""
    # u_ft0 for default Dp=75µm, rho_a=1.2 ≈ 0.206 m/s; set ustar well below
    c = {**cell, "ustar": np.array([[0.05]])}
    assert dust_emission(**c)[0, 0] == 0.0


def test_zero_flux_full_vegetation(cell):
    """Emission must be zero when VAI >= vai_thr (f_bare = 0)."""
    c = {**cell, "vai": np.array([[VAI_THR]])}
    assert dust_emission(**c)[0, 0] == 0.0


def test_zero_flux_zero_clay(cell):
    """Emission must be zero when f_clay = 0 (fclay0 = 0, no sandblasting)."""
    c = {**cell, "f_clay": np.array([[0.0]])}
    assert dust_emission(**c)[0, 0] == 0.0


# --- absolute value / dimensional check ---


def test_flux_hand_calc(cell):
    """F_d must match an independent hand calculation of the K14a equation
    (guards against dimensional errors that relative tests cannot catch)."""
    rho_a, ustar = 1.2, 0.5
    u_t = np.sqrt(
        SL_A * (RHOP_DEFAULT * GRAV * DP_DEFAULT + SL_GAMMA / DP_DEFAULT) / rho_a
    )  # f_m = 1 (dry)
    u_st = u_t * np.sqrt(rho_a / RHO_A0)
    C_d = CD0 * np.exp(-CE * (u_st - U_ST0) / U_ST0)
    kappa = min(C_KAPPA * (u_st - U_ST0) / U_ST0, KAPPA_MAX)
    f_bare = 1.0 - 0.05 / VAI_THR
    expected = (
        C_TUNE_DEFAULT
        * C_d
        * f_bare
        * 0.15
        * rho_a
        * (ustar**2 - u_t**2)
        / u_st
        * (ustar / u_t) ** kappa
    )
    assert np.isclose(dust_emission(**cell)[0, 0], expected, rtol=1e-12)
    # order-of-magnitude anchor for the demo forcing
    assert np.isclose(expected, 3.737e-7, rtol=2e-3)


# --- monotonicity and scaling ---


def test_flux_increases_with_ustar(cell):
    """Higher ustar above threshold must produce larger flux."""
    F_low = dust_emission(**{**cell, "ustar": np.array([[0.4]])})
    F_high = dust_emission(**{**cell, "ustar": np.array([[0.8]])})
    assert F_high[0, 0] > F_low[0, 0]


def test_flux_decreases_with_soil_moisture(cell):
    """Higher soil moisture raises u_ft, reducing or zeroing emission."""
    F_dry = dust_emission(**{**cell, "w": np.array([[0.01]])})
    F_moist = dust_emission(**{**cell, "w": np.array([[0.15]])})
    assert F_moist[0, 0] <= F_dry[0, 0]


def test_flux_scales_linearly_with_C_tune(cell):
    """Doubling C_tune must exactly double the emission flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, "C_tune": cell["C_tune"] * 2.0})
    assert np.isclose(F2[0, 0], F1[0, 0] * 2.0, rtol=1e-10)


def test_flux_increases_with_f_clay_below_cap(cell):
    """Flux must increase with f_clay when f_clay < FCLAY_MAX."""
    F_low = dust_emission(**{**cell, "f_clay": np.array([[0.05]])})
    F_high = dust_emission(**{**cell, "f_clay": np.array([[0.15]])})
    assert F_high[0, 0] > F_low[0, 0]


def test_flux_flat_above_clay_cap(cell):
    """Flux must be identical when f_clay is above FCLAY_MAX (cap enforced)."""
    F1 = dust_emission(**{**cell, "f_clay": np.array([[0.25]])})
    F2 = dust_emission(**{**cell, "f_clay": np.array([[0.40]])})
    assert np.isclose(F1[0, 0], F2[0, 0], rtol=1e-10)


# --- Leung et al. (2023) parameter options ---


def test_leung2023_Dp_lowers_threshold(cell):
    """Using Dp=127µm (Leung 2023) vs Dp=75µm changes u_ft0."""
    u_k14 = threshold_velocity_dry(Dp=75e-6, rho_a=cell["rho_a"])
    u_leung = threshold_velocity_dry(Dp=127e-6, rho_a=cell["rho_a"])
    # Shao & Lu parabola: minimum near 80µm; 127µm is on the larger-particle
    # side, so u_ft0(127µm) > u_ft0(75µm)
    assert u_leung[0, 0] > u_k14[0, 0]


def test_higher_vai_thr_increases_f_bare(cell):
    """Using vai_thr=1.0 (Leung 2023) gives larger f_bare for moderate VAI."""
    vai = np.array([[0.15]])
    f_k14 = bare_fraction(vai, vai_thr=0.3)
    f_leung = bare_fraction(vai, vai_thr=1.0)
    assert f_leung[0, 0] > f_k14[0, 0]


# --- vectorized grid ---


def test_vectorized_grid():
    """Must run correctly on a full 2D grid with mixed land/ocean cells."""
    nx, ny = 12, 10
    rng = np.random.default_rng(42)
    F_d = dust_emission(
        ustar=rng.uniform(0.1, 1.0, (nx, ny)),
        w=rng.uniform(0.0, 0.3, (nx, ny)),
        rho_a=rng.uniform(1.0, 1.3, (nx, ny)),
        vai=rng.uniform(0.0, 0.5, (nx, ny)),
        snowdp=rng.uniform(0.0, 0.02, (nx, ny)),
        f_clay=rng.uniform(0.05, 0.4, (nx, ny)),
        oro=rng.choice([LAND, 2.0], (nx, ny)),
    )
    assert F_d.shape == (nx, ny)
    assert (F_d >= 0.0).all()


# --- K14 constant guards ---


def test_cd0_value():
    """CD0 must match Kok et al. (2014a) Table 1: (4.4 ± 0.5) × 10^{-5}."""
    assert np.isclose(CD0, 4.4e-5, rtol=0.01), f"CD0 = {CD0:.2e}; expected 4.4e-5"


def test_ce_value():
    """CE must match Kok et al. (2014a) Table 1: 2.0 ± 0.3."""
    assert np.isclose(CE, 2.0, atol=0.01), f"CE = {CE}; expected 2.0"


def test_u_st0_value():
    """U_ST0 must match Kok et al. (2012): 0.16 m/s."""
    assert np.isclose(U_ST0, 0.16, atol=1e-4), f"U_ST0 = {U_ST0}; expected 0.16"


def test_c_kappa_value():
    """C_KAPPA must match Kok et al. (2014a): 2.7 ± 1.0."""
    assert np.isclose(C_KAPPA, 2.7, atol=0.01), f"C_KAPPA = {C_KAPPA}; expected 2.7"


def test_kappa_max_value():
    """KAPPA_MAX must equal 3.0 per Leung et al. (2023)."""
    assert np.isclose(KAPPA_MAX, 3.0, atol=1e-10), f"KAPPA_MAX = {KAPPA_MAX}; expected 3.0"


# --- initialize() ---


def test_initialize_missing_file():
    """initialize() must raise FileNotFoundError for missing paths."""
    with pytest.raises(FileNotFoundError):
        initialize("no_fclay.npy", "no_oro.npy")


def test_initialize_fclay_wrong_range(tmp_path):
    """initialize() must raise ValueError if f_clay values exceed [0, 1]."""
    fc_f = tmp_path / "fclay.npy"
    or_f = tmp_path / "oro.npy"
    np.save(fc_f, np.ones((4, 4)) * 1.5)  # out of range
    np.save(or_f, np.ones((4, 4)))
    with pytest.raises(ValueError, match="f_clay"):
        initialize(str(fc_f), str(or_f))


def test_initialize_shape_mismatch(tmp_path):
    """initialize() must raise ValueError if f_clay and oro shapes differ."""
    fc_f = tmp_path / "fclay.npy"
    or_f = tmp_path / "oro.npy"
    np.save(fc_f, np.ones((4, 4)) * 0.15)
    np.save(or_f, np.ones((4, 6)))  # wrong shape
    with pytest.raises(ValueError, match="oro shape"):
        initialize(str(fc_f), str(or_f))


def test_initialize_valid(tmp_path):
    """initialize() must return correct shapes for valid inputs."""
    fc_f = tmp_path / "fclay.npy"
    or_f = tmp_path / "oro.npy"
    np.save(fc_f, np.ones((6, 8)) * 0.15)
    np.save(or_f, np.ones((6, 8)))
    f_clay, oro = initialize(str(fc_f), str(or_f))
    assert f_clay.shape == (6, 8)
    assert oro.shape == (6, 8)
