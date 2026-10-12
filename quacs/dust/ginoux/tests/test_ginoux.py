"""Tests for the GOCART dust emission scheme (Ginoux et al., 2001)."""

import numpy as np
import pytest

from quacs.dust.ginoux.ginoux import (
    C_DEFAULT,
    LAND,
    RADIUS_DEFAULT,
    RHOP_DEFAULT,
    SP_DEFAULT,
    dust_emission,
    initialize,
    threshold_velocity,
)

# ---------------------------------------------------------------------------
# Pytest test suite
# Run with: pytest quacs/dust/ginoux -v
# ---------------------------------------------------------------------------


@pytest.fixture
def cell():
    """Single land cell, moderate moisture, wind above threshold."""
    return dict(
        w10m=np.array([[12.0]]),
        u_t=np.array([[6.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]),
        S=np.array([[0.5]]),
        s_p=SP_DEFAULT.copy(),
        C=C_DEFAULT,
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
        ratio_flux = F_p[i] / F_p[i + 1]
        ratio_sp = cell["s_p"][i] / cell["s_p"][i + 1]
        assert np.isclose(ratio_flux, ratio_sp, rtol=1e-10), (
            f"Bin {i}/{i + 1} flux ratio {ratio_flux:.4f} != s_p ratio {ratio_sp:.4f}"
        )


def test_sp_default_matches_fortran():
    """
    SP_DEFAULT must match DU2G_instance.F90, L43 exactly.
    Note: source_fraction sums to 1.1 in the Fortran; this is expected.
    See module docstring for explanation of why s_p is not a strict partition.
    """
    expected = np.array([0.1, 0.25, 0.25, 0.25, 0.25])
    assert np.allclose(SP_DEFAULT, expected), (
        f"SP_DEFAULT {SP_DEFAULT} does not match Fortran source_fraction {expected}"
    )
    assert np.isclose(SP_DEFAULT.sum(), 1.1), (
        f"Expected sum 1.1 (as in Fortran); got {SP_DEFAULT.sum():.4f}"
    )


# --- analytic correctness ---


def test_total_flux_matches_analytic(cell):
    """Total flux must match the analytic formula exactly."""
    c = cell
    gwet_clamped = max(c["gwet"][0, 0], 1e-3)
    u_thresh = max(c["u_t"][0, 0] * (1.2 + 0.2 * np.log10(gwet_clamped)), 0.0)
    # base_flux * sum(s_p) = total across bins
    base_flux = c["C"] * c["S"][0, 0] * c["w10m"][0, 0] ** 2 * (c["w10m"][0, 0] - u_thresh)
    expected = base_flux * c["s_p"].sum()
    F_p = dust_emission(**c)
    assert np.isclose(F_p[0, 0, :].sum(), expected, rtol=1e-10)


# --- emission cutoffs: land mask ---


def test_zero_flux_over_ocean():
    """Emission must be zero over non-land cells."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]),
        u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[2.0]]),
        frlake=np.array([[0.0]]),
        S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()


# --- emission cutoffs: soil moisture ---


def test_zero_flux_wet_soil():
    """Emission must be zero when gwet >= 0.5."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]),
        u_t=np.array([[5.0]]),
        gwet=np.array([[0.5]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]),
        S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()


def test_moisture_correction_raises_threshold():
    """Moister soil raises effective threshold, reducing flux."""
    base = dict(
        w10m=np.array([[10.0]]),
        u_t=np.array([[6.0]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]),
        S=np.array([[0.5]]),
        s_p=SP_DEFAULT,
    )
    F_dry = dust_emission(gwet=np.array([[0.01]]), **base)
    F_moist = dust_emission(gwet=np.array([[0.40]]), **base)
    assert F_moist[0, 0, :].sum() < F_dry[0, 0, :].sum()


# --- emission cutoffs: wind threshold ---


def test_zero_flux_below_threshold():
    """Wind below adjusted threshold must give zero flux."""
    F_p = dust_emission(
        w10m=np.array([[3.0]]),
        u_t=np.array([[8.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]),
        S=np.array([[1.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()


# --- lake fraction ---


def test_full_lake_gives_zero_flux():
    """frlake = 1.0 must give zero flux."""
    F_p = dust_emission(
        w10m=np.array([[12.0]]),
        u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[1.0]]),
        S=np.array([[0.5]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()


def test_lake_fraction_scales_flux():
    """frlake = 0.5 must give exactly half the flux of frlake = 0."""
    base = dict(
        w10m=np.array([[12.0]]),
        u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        S=np.array([[0.5]]),
        s_p=SP_DEFAULT,
    )
    F_no_lake = dust_emission(frlake=np.array([[0.0]]), **base)
    F_half_lake = dust_emission(frlake=np.array([[0.5]]), **base)
    assert np.isclose(F_half_lake.sum(), F_no_lake.sum() * 0.5, rtol=1e-10)


# --- zero erodibility ---


def test_zero_erodibility_gives_zero_flux():
    """S = 0 must give zero flux regardless of all other conditions."""
    F_p = dust_emission(
        w10m=np.array([[15.0]]),
        u_t=np.array([[5.0]]),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        frlake=np.array([[0.0]]),
        S=np.array([[0.0]]),
        s_p=SP_DEFAULT,
    )
    assert (F_p == 0.0).all()


# --- scaling laws ---


def test_flux_scales_linearly_with_C(cell):
    """Doubling C must exactly double the total flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, "C": cell["C"] * 2.0})
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-10)


def test_flux_scales_linearly_with_S(cell):
    """Doubling S must exactly double the total flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, "S": cell["S"] * 2.0})
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-10)


# --- vectorized grid ---


def test_vectorized_grid():
    """Must run correctly on a full 2D grid with mixed land/ocean cells."""
    nx, ny, nb = 12, 10, 5
    rng = np.random.default_rng(0)
    F_p = dust_emission(
        w10m=rng.uniform(0, 20, (nx, ny)),
        u_t=rng.uniform(4, 8, (nx, ny)),
        gwet=rng.uniform(0, 0.45, (nx, ny)),
        oro=rng.choice([LAND, 2.0], (nx, ny)),
        frlake=rng.uniform(0, 1, (nx, ny)),
        S=rng.uniform(0, 1, (nx, ny)),
        s_p=SP_DEFAULT,
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
    assert np.all(np.diff(u_t[0, 0, :]) < 0), (
        f"Expected decreasing threshold across bins; got {u_t[0, 0, :]}"
    )


# --- source-verified parameter guards ---


def test_radius_bin1_is_073_micron():
    """
    Bin 1 radius must be 0.73 um (DU2G_instance.F90, L26).
    Guards against reverting to the earlier placeholder value of 0.5 um.
    """
    assert np.isclose(RADIUS_DEFAULT[0], 0.73e-6), (
        f"Bin 1 radius should be 0.73e-6 m; got {RADIUS_DEFAULT[0]:.2e}"
    )


def test_particle_density_per_bin():
    """
    Bin 1 density must be 2500 kg/m^3; bins 2-5 must be 2650 kg/m^3.
    Source: DU2G_instance.F90, L33.
    """
    assert RHOP_DEFAULT[0] == 2500.0, f"Bin 1 density: expected 2500, got {RHOP_DEFAULT[0]}"
    assert np.all(RHOP_DEFAULT[1:] == 2650.0), (
        f"Bins 2-5 density: expected 2650, got {RHOP_DEFAULT[1:]}"
    )


# --- initialize() ---


def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("no_S.npy", "no_sp.npy", "no_oro.npy", "no_frlake.npy")


def test_initialize_sp_wrong_shape(tmp_path):
    """initialize() must raise ValueError if s_p has wrong shape."""
    S_f = tmp_path / "S.npy"
    sp_f = tmp_path / "sp.npy"
    or_f = tmp_path / "oro.npy"
    fl_f = tmp_path / "frl.npy"
    np.save(S_f, np.ones((4, 4)) * 0.5)
    np.save(sp_f, np.array([0.25, 0.25, 0.25, 0.25]))  # 4 bins, expected 5
    np.save(or_f, np.ones((4, 4)))
    np.save(fl_f, np.zeros((4, 4)))
    with pytest.raises(ValueError, match="shape"):
        initialize(str(S_f), str(sp_f), str(or_f), str(fl_f))


def test_initialize_valid(tmp_path):
    """initialize() must return correct shapes for valid inputs."""
    S_f = tmp_path / "S.npy"
    sp_f = tmp_path / "sp.npy"
    or_f = tmp_path / "oro.npy"
    fl_f = tmp_path / "frl.npy"
    np.save(S_f, np.ones((6, 6)) * 0.4)
    np.save(sp_f, SP_DEFAULT)
    np.save(or_f, np.ones((6, 6)))
    np.save(fl_f, np.zeros((6, 6)))
    S, s_p, oro, frlake = initialize(str(S_f), str(sp_f), str(or_f), str(fl_f))
    assert S.shape == (6, 6)
    assert s_p.shape == (5,)
    assert oro.shape == (6, 6)
    assert frlake.shape == (6, 6)
