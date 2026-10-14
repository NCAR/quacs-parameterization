"""Tests for quacs.dust.afwa."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import quacs.dust.afwa.afwa as afwa
from quacs.dust.afwa.afwa import (
    ALPHA_DEFAULT,
    BETAMAX,
    CMB,
    DEN_SALT,
    DISTR_DUST,
    FRAC_SALT,
    GAMMA_DEFAULT,
    GRAV,
    LAND,
    NDUST,
    NGSALT,
    NONERODIBLE_SOIL,
    POROSITY_TABLE,
    REFF_SALT,
    SPOINT,
    AfwaSettings,
    afwa_dust_driver,
    dust_emission_afwa,
    initialize,
    saltation_area_fraction,
    threshold_friction_velocity,
    tunable_names,
)


def _afwa_cell():
    """Single erodible land cell, dry soil, strong friction velocity."""
    return dict(
        ustar=np.array([[0.6]]),
        volsm=np.array([[0.02]]),
        airden=np.array([[1.225]]),
        znt=np.array([[0.01]]),
        snowh=np.array([[0.0]]),
        erod=np.array([[[0.8, 0.0, 0.0]]]),
        clay=np.array([[0.20]]),
        sand=np.array([[0.60]]),
        oro=np.array([[LAND]]),
        isltyp=np.array([[3]]),
    )


@pytest.fixture
def cell():
    return _afwa_cell()


# --- shape and basic contract ---


def test_output_shape(cell):
    """F_p must have shape [x, y, n_dust]."""
    assert dust_emission_afwa(**cell).shape == (1, 1, NDUST)


def test_all_fluxes_nonnegative(cell):
    """No negative emission fluxes anywhere."""
    assert (dust_emission_afwa(**cell) >= 0.0).all()


def test_invalid_land_texture_warns_but_preserves_wrf_result(cell):
    """Unphysical land texture needs a warning without changing WRF arithmetic."""
    altered = {**cell, "clay": np.array([[0.6]]), "sand": np.array([[0.6]])}
    with pytest.warns(UserWarning, match="clay|sand"):
        result = dust_emission_afwa(**altered)
    assert result.shape == (1, 1, NDUST)


def test_invalid_sea_texture_does_not_warn(cell):
    """WRF skips water cells before reading their texture inputs."""
    import warnings

    altered = {
        **cell,
        "oro": np.array([[2.0]]),
        "clay": np.array([[-1.0]]),
        "sand": np.array([[2.0]]),
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        assert (dust_emission_afwa(**altered) == 0).all()


def test_invalid_land_erod_warns_without_rejecting(cell):
    """Out-of-range land erodibility is reported while WRF still computes."""
    altered = {**cell, "erod": np.array([[[1.2, 0.0, 0.0]]])}
    with pytest.warns(UserWarning, match="erod"):
        assert dust_emission_afwa(**altered).sum() > 0


def test_invalid_land_volsm_warns_without_rejecting(cell):
    """Volumetric moisture above one is diagnosed but follows WRF physics."""
    altered = {**cell, "volsm": np.array([[1.2]])}
    with pytest.warns(UserWarning, match="volsm"):
        assert dust_emission_afwa(**altered).shape == (1, 1, NDUST)


def test_missing_selected_dataset_stays_an_error_over_sea(cell):
    """A selected regional dataset remains required for an all-sea grid."""
    altered = {**cell, "oro": np.array([[2.0]])}
    with pytest.raises(ValueError, match="erod_dri"):
        dust_emission_afwa(**altered, dust_dsr=1)


def test_wrong_selected_dataset_shape_stays_an_error(cell):
    """A chosen regional erodibility array must match the grid and bins."""
    with pytest.raises(ValueError, match="erod_dri"):
        dust_emission_afwa(**cell, dust_dsr=1, erod_dri=np.ones((2, 1, 3)))


def test_emission_is_positive_for_strong_wind(cell):
    """A dry, erodible, windy cell must emit."""
    assert dust_emission_afwa(**cell).sum() > 0.0


# --- bin distribution ---


def test_bin_distribution_matches_distr_dust(cell):
    """Each dust bin must carry exactly distr_dust[m] of the total emission."""
    F_p = dust_emission_afwa(**cell)[0, 0, :]
    total = F_p.sum()
    assert np.allclose(F_p / total, DISTR_DUST, rtol=1e-10)


def test_distr_dust_sums_to_one():
    """distr_dust (Kok 2011) must sum to 1.0 (source: module_data_gocart_dust.F)."""
    assert np.isclose(DISTR_DUST.sum(), 1.0, atol=1e-4)


def test_total_flux_equals_emit(cell):
    """Sum over dust bins equals the total emitted flux (distr_dust sums to 1)."""
    F_p = dust_emission_afwa(**cell)
    # recompute emit independently
    us = cell["ustar"][0, 0]
    rhoa = cell["airden"][0, 0]
    u_ts0 = threshold_friction_velocity(REFF_SALT, DEN_SALT, cell["airden"])[0, 0, :]
    ds_rel = saltation_area_fraction(cell["clay"], cell["sand"])[0, 0, :]
    salt = np.where(
        us > u_ts0,
        CMB * ds_rel * (rhoa / GRAV) * us**3 * (1.0 + u_ts0 / us) * (1.0 - u_ts0**2 / us**2),
        0.0,
    )
    beta = min(10.0 ** (13.6 * cell["clay"][0, 0] - 6.0), BETAMAX)
    emit = salt.sum() * cell["erod"][0, 0].sum() ** GAMMA_DEFAULT * ALPHA_DEFAULT * beta
    assert np.isclose(F_p.sum(), emit, rtol=1e-9)


# --- emission cutoffs ---


def test_zero_flux_over_ocean(cell):
    """Emission must be zero over non-land cells."""
    F_p = dust_emission_afwa(**{**cell, "oro": np.array([[2.0]])})
    assert (F_p == 0.0).all()


def test_zero_flux_zero_erodibility(cell):
    """erod = 0 must give zero flux."""
    F_p = dust_emission_afwa(**{**cell, "erod": np.array([[[0.0, 0.0, 0.0]]])})
    assert (F_p == 0.0).all()


def test_zero_flux_below_threshold(cell):
    """u* below the threshold friction velocity must give zero flux."""
    F_p = dust_emission_afwa(**{**cell, "ustar": np.array([[0.05]])})
    assert (F_p == 0.0).all()


def test_zero_flux_rough_surface(cell):
    """Roughness length above ZNT_MAX must suppress emission."""
    F_p = dust_emission_afwa(**{**cell, "znt": np.array([[0.5]])})
    assert (F_p == 0.0).all()


def test_zero_flux_snow_cover(cell):
    """Snow depth above SNOW_MAX must suppress emission."""
    F_p = dust_emission_afwa(**{**cell, "snowh": np.array([[0.05]])})
    assert (F_p == 0.0).all()


def test_zero_flux_vegetated(cell):
    """A zero LAI mask (fully vegetated) must give zero flux."""
    F_p = dust_emission_afwa(**{**cell, "dust_veg": 2, "lai_vegmask": np.array([[0.0]])})
    assert (F_p == 0.0).all()


def test_zero_flux_nonerodible_soil(cell):
    """Bedrock / land-ice / lava soil types must give zero flux."""
    for st in NONERODIBLE_SOIL:
        F_p = dust_emission_afwa(**{**cell, "isltyp": np.array([[st]])})
        assert (F_p == 0.0).all(), f"soil type {st} should not loft"


# --- soil moisture ---


def test_moisture_raises_threshold_reduces_flux(cell):
    """Wetter soil (above dry limit) must raise the threshold and reduce flux."""
    F_dry = dust_emission_afwa(**{**cell, "volsm": np.array([[0.02]])}).sum()
    F_wet = dust_emission_afwa(**{**cell, "volsm": np.array([[0.25]])}).sum()
    assert F_wet < F_dry


# --- scaling laws ---


def test_flux_scales_linearly_with_alpha(cell):
    """Doubling alpha must exactly double the total flux."""
    F1 = dust_emission_afwa(**cell).sum()
    F2 = dust_emission_afwa(**cell, settings=AfwaSettings(alpha=2.0 * ALPHA_DEFAULT)).sum()
    assert np.isclose(F2, F1 * 2.0, rtol=1e-10)


def test_higher_erodibility_increases_flux(cell):
    """With gamma = 1, higher erodibility must increase flux monotonically."""
    F_lo = dust_emission_afwa(**{**cell, "erod": np.array([[[0.3, 0.0, 0.0]]])}).sum()
    F_hi = dust_emission_afwa(**{**cell, "erod": np.array([[[0.9, 0.0, 0.0]]])}).sum()
    assert F_hi > F_lo


# --- sandblasting efficiency beta ---


def test_beta_formula_below_cap():
    """beta = 10^(13.6*clay - 6) below the cap (clay = 0.05)."""
    clay = 0.05
    expected = 10.0 ** (13.6 * clay - 6.0)
    assert expected < BETAMAX
    beta = min(10.0 ** (13.6 * clay - 6.0), BETAMAX)
    assert np.isclose(beta, expected)


def test_beta_saturates_at_betamax():
    """High clay (>= ~0.20) must saturate beta at BETAMAX (afwa.F)."""
    clay = 0.30
    beta = min(10.0 ** (13.6 * clay - 6.0), BETAMAX)
    assert np.isclose(beta, BETAMAX)


# --- threshold_friction_velocity ---


def test_threshold_shape():
    """threshold_friction_velocity must return shape [x, y, smx]."""
    u = threshold_friction_velocity(REFF_SALT, DEN_SALT, np.ones((4, 6)) * 1.2)
    assert u.shape == (4, 6, NGSALT)


def test_threshold_positive():
    """All threshold friction velocities must be positive."""
    u = threshold_friction_velocity(REFF_SALT, DEN_SALT, np.ones((3, 3)) * 1.2)
    assert (u > 0.0).all()


def test_threshold_is_u_shaped():
    """
    Across the 9 saltation bins (0.71-125 um) the threshold is U-shaped:
    it has an interior minimum (Iversen-White), not a monotonic trend.
    """
    u = threshold_friction_velocity(REFF_SALT, DEN_SALT, np.ones((1, 1)) * 1.225)[0, 0, :]
    imin = int(np.argmin(u))
    assert 0 < imin < NGSALT - 1, f"expected interior minimum, got bin {imin}"


# --- saltation_area_fraction ---


def test_ds_rel_sums_to_one():
    """Relative surface-area fractions must sum to 1 over saltation bins."""
    ds = saltation_area_fraction(np.array([[0.2]]), np.array([[0.6]]))
    assert np.isclose(ds.sum(), 1.0, rtol=1e-12)


# --- source-verified parameter guards ---


def test_ngsalt_is_nine():
    """ngsalt must be 9 (module_data_gocart_dust.F)."""
    assert NGSALT == 9
    assert REFF_SALT.shape == (9,) and DEN_SALT.shape == (9,)
    assert SPOINT.shape == (9,) and FRAC_SALT.shape == (9,)


def test_spoint_matches_source():
    """spoint must be (1,2,2,2,2,2,3,3,3) (module_data_gocart_dust.F)."""
    assert np.array_equal(SPOINT, np.array([1, 2, 2, 2, 2, 2, 3, 3, 3]))


def test_reff_salt_endpoints():
    """First/last saltation radii must be 0.71 um and 125 um."""
    assert np.isclose(REFF_SALT[0], 0.71e-6)
    assert np.isclose(REFF_SALT[-1], 125.0e-6)


def test_distr_dust_matches_source():
    """distr_dust must match the hardcoded Kok values in module_data_gocart_dust.F."""
    expected = np.array([1.074e-1, 1.012e-1, 2.078e-1, 4.817e-1, 1.019e-1])
    assert np.allclose(DISTR_DUST, expected)


def test_betamax_and_cmb():
    """betamax and active cmb must match afwa.F."""
    assert np.isclose(BETAMAX, 5.25e-4)
    assert np.isclose(CMB, 1.0)  # White (1979) 2.61 is commented out in source


# --- vectorized grid ---


def test_vectorized_grid():
    """Must run correctly on a full 2D grid with mixed land/ocean cells."""
    nx, ny = 12, 10
    rng = np.random.default_rng(0)
    F_p = dust_emission_afwa(
        ustar=rng.uniform(0.0, 1.0, (nx, ny)),
        volsm=rng.uniform(0.0, 0.3, (nx, ny)),
        airden=rng.uniform(1.0, 1.3, (nx, ny)),
        znt=rng.uniform(0.0, 0.3, (nx, ny)),
        snowh=rng.uniform(0.0, 0.05, (nx, ny)),
        erod=rng.uniform(0.0, 1.0, (nx, ny, 3)),
        clay=rng.uniform(0.0, 0.4, (nx, ny)),
        sand=rng.uniform(0.3, 0.7, (nx, ny)),
        oro=rng.choice([LAND, 2.0], (nx, ny)),
        isltyp=np.full((nx, ny), 3),
    )
    assert F_p.shape == (nx, ny, NDUST)
    assert (F_p >= 0.0).all()
    assert np.isfinite(F_p).all()


# --- initialize() ---


def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("no_e.npy", "no_c.npy", "no_s.npy", "no_i.npy", "no_o.npy")


def test_initialize_valid(tmp_path):
    """initialize() retains WRF soil types for the timestep porosity lookup."""
    paths = {n: tmp_path / f"{n}.npy" for n in ["erod", "clay", "sand", "isltyp", "oro"]}
    np.save(paths["erod"], np.ones((6, 6, 3)) * 0.5)
    np.save(paths["clay"], np.ones((6, 6)) * 0.2)
    np.save(paths["sand"], np.ones((6, 6)) * 0.6)
    np.save(paths["isltyp"], np.full((6, 6), 3))  # -> porosity 0.434
    np.save(paths["oro"], np.ones((6, 6)))
    erod, clay, sand, soil, oro, vegfra, lai, dri, cnga, snga = initialize(
        str(paths["erod"]),
        str(paths["clay"]),
        str(paths["sand"]),
        str(paths["isltyp"]),
        str(paths["oro"]),
    )
    assert erod.shape == (6, 6, 3)
    assert np.allclose(soil, 3)  # WRF soil type 3 is retained
    assert vegfra is None and lai is None


def test_initialize_preserves_clay_plus_sand_above_one(tmp_path):
    """WRF texture inputs must be loaded even when clay + sand exceeds one."""
    paths = {n: tmp_path / f"{n}.npy" for n in ["erod", "clay", "sand", "isltyp", "oro"]}
    np.save(paths["erod"], np.ones((4, 4, 3)) * 0.5)
    np.save(paths["clay"], np.ones((4, 4)) * 0.6)
    np.save(paths["sand"], np.ones((4, 4)) * 0.6)  # 0.6 + 0.6 > 1
    np.save(paths["isltyp"], np.full((4, 4), 3))
    np.save(paths["oro"], np.ones((4, 4)))
    result = initialize(
        str(paths["erod"]),
        str(paths["clay"]),
        str(paths["sand"]),
        str(paths["isltyp"]),
        str(paths["oro"]),
    )
    assert result[1][0, 0] + result[2][0, 0] == 1.2


# --- AFWA regression tests (A1-A5, A9, A11, A12) ---


def test_port_gravity_preserves_wrf_flux(cell):
    """A1: WRF gravity must reproduce the pinned 8-byte emission without overrides."""
    expected = np.array(
        [
            6.58893048854812215e-8,
            6.20856392403230911e-8,
            1.27484148558687142e-7,
            2.95520280850431189e-7,
            6.25150853615506318e-8,
        ]
    )
    np.testing.assert_allclose(dust_emission_afwa(**cell)[0, 0], expected, rtol=1e-12, atol=0)


def test_port_fractional_land_preserves_wrf_emission(cell):
    """A2: Fractional xland below 1.5 must retain land emission."""
    reference = dust_emission_afwa(**cell)
    for land in [1.25, 1.49]:
        np.testing.assert_array_equal(
            dust_emission_afwa(**{**cell, "oro": np.array([[land]])}), reference
        )
    for land in [1.5, 2.0]:
        assert (dust_emission_afwa(**{**cell, "oro": np.array([[land]])}) == 0).all()


def test_port_initializer_documents_unmodified_erodibility():
    """A3: Callers must know initialization returns raw erodibility."""
    doc = initialize.__doc__
    assert "raw erodibility" in doc and "does not apply" in doc


def test_port_required_soil_mask_documents_wrf_exclusions():
    """A4: Callers must know the required soil mask preserves WRF exclusions."""
    doc = dust_emission_afwa.__doc__
    assert "15, 16, 18" in doc and "isltyp is required" in doc
    assert "238-241" in doc


def test_port_required_vegetation_documents_ginoux_default():
    """A5: Callers must know an omitted mask cannot recover Ginoux bin 1."""
    doc = dust_emission_afwa.__doc__
    assert "Ginoux" in doc and "internally constructed mask" in doc
    assert "185-191" in doc


def test_port_invalid_area_does_not_silently_suppress_dust(cell):
    """A9: Direct Fortran division must not turn invalid negative area into zero."""
    invalid = {
        **cell,
        "clay": np.array([[-1.0]]),
        "sand": np.array([[0.0]]),
        "volsm": np.array([[0.0]]),
        "ustar": np.array([[10.0]]),
    }
    fractions = saltation_area_fraction(invalid["clay"], invalid["sand"])
    assert np.isclose(fractions.sum(), 1.0)
    assert (dust_emission_afwa(**invalid) > 0).all()


def test_port_calm_wind_is_warning_free(cell):
    """A11: Calm cells must skip saltation divisions and keep emitting cells exact."""
    mixed = {k: np.repeat(v, 2, axis=1) for k, v in cell.items()}
    mixed["ustar"][0, 0] = 0.0
    with np.errstate(divide="raise", invalid="raise"):
        flux = dust_emission_afwa(**mixed)
    assert (flux[:, 0] == 0).all()
    np.testing.assert_array_equal(flux[:, 1], dust_emission_afwa(**cell)[:, 0])


def test_port_numpy_only_import_keeps_scheme_usable():
    """A12: Model users must be able to import and run without pytest installed."""
    import subprocess
    import sys
    from textwrap import dedent

    code = """
import builtins
import runpy
import sys
import numpy as np
original = builtins.__import__
def without_pytest(name, *args, **kwargs):
    if name == 'pytest' or name.startswith('pytest.'):
        raise ImportError('pytest unavailable in model environment')
    return original(name, *args, **kwargs)
builtins.__import__ = without_pytest
scheme = runpy.run_path(sys.argv[1], run_name='afwa_dependency_probe')
a = lambda x: np.array([[x]])
flux = scheme['dust_emission_afwa'](
    a(.6), a(.02), a(1.225), a(.01), a(0),
    np.array([[[.8, 0, 0]]]), a(.2), a(.6), a(1), a(3))
assert flux.shape == (1, 1, 5) and (flux > 0).all()
"""
    result = subprocess.run(
        [sys.executable, "-B", "-c", dedent(code), str(Path(afwa.__file__))],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


# --- WRF soil and vegetation gates ---


def test_a4a5_missing_soil_type_cannot_skip_wrf_exclusions(cell):
    """A4: Missing isltyp must fail before excluded soils can emit."""
    fields = dict(cell)
    fields.pop("isltyp", None)
    with pytest.raises(TypeError, match="isltyp"):
        dust_emission_afwa(**fields)


def test_a4a5_bedrock_ice_and_lava_cannot_emit(cell):
    """A4: WRF soil types 15, 16 and 18 must suppress wind-driven dust."""
    for soil in [15, 16, 18]:
        fields = {**cell, "erod": np.array([[[0.8, 0.0, 0.0]]]), "isltyp": np.array([[soil]])}
        assert (dust_emission_afwa(**fields) == 0).all()


def test_a4a5_ginoux_first_bin_controls_default_vegetation(cell):
    """A5: Nonzero total erodibility cannot override a zero Ginoux first bin."""
    fields = {
        **cell,
        "erod": np.array([[[0.0, 0.3, 0.5]]]),
        "isltyp": np.array([[3]]),
        "dust_veg": 0,
    }
    assert (dust_emission_afwa(**fields) == 0).all()


def test_a4a5_greenfrac_five_percent_suppresses_dust(cell):
    """A5: The WRF green-fraction cutoff is inclusive at five percent."""
    fields = {
        **cell,
        "erod": np.array([[[0.8, 0.0, 0.0]]]),
        "isltyp": np.array([[3]]),
        "dust_veg": 1,
    }
    assert (dust_emission_afwa(**{**fields, "vegfra": np.array([[5.0]])}) == 0).all()
    assert (dust_emission_afwa(**{**fields, "vegfra": np.array([[4.9]])}) > 0).all()


def test_a4a5_lai_mask_suppresses_vegetated_dust(cell):
    """A5: A zero MODIS LAI mask must suppress otherwise erodible soil."""
    fields = {
        **cell,
        "erod": np.array([[[0.8, 0.0, 0.0]]]),
        "isltyp": np.array([[3]]),
        "dust_veg": 2,
        "lai_vegmask": np.array([[0.0]]),
    }
    assert (dust_emission_afwa(**fields) == 0).all()


def test_a4a5_greenfrac_option_requires_vegetation_data(cell):
    """A5: Selecting green fraction must not silently invent its dataset."""
    with pytest.raises(ValueError, match="vegfra"):
        dust_emission_afwa(
            **{
                **cell,
                "erod": np.array([[[0.8, 0.0, 0.0]]]),
                "isltyp": np.array([[3]]),
                "dust_veg": 1,
            }
        )


def test_a4a5_lai_option_requires_vegetation_data(cell):
    """A5: Selecting LAI must not silently invent its mask dataset."""
    with pytest.raises(ValueError, match="lai_vegmask"):
        dust_emission_afwa(
            **{
                **cell,
                "erod": np.array([[[0.8, 0.0, 0.0]]]),
                "isltyp": np.array([[3]]),
                "dust_veg": 2,
            }
        )


def test_a4a5_unknown_option_retains_wrf_default_mask(cell):
    """New 2: An unknown selector must preserve the WRF Ginoux fallback."""
    np.testing.assert_array_equal(
        dust_emission_afwa(**cell, dust_veg=3), dust_emission_afwa(**cell)
    )


def test_a4a5_initialization_loads_only_selected_vegetation(tmp_path):
    """A5: Initialization must load the selected WRF dataset without inventing a mask."""
    paths = [tmp_path / (name + ".npy") for name in ["erod", "clay", "sand", "isltyp", "oro"]]
    values = [
        np.array([[[0.8, 0.0, 0.0]]]),
        np.array([[0.2]]),
        np.array([[0.6]]),
        np.array([[3]]),
        np.array([[1.0]]),
    ]
    for path, value in zip(paths, values):
        np.save(path, value)
    green = tmp_path / "green.npy"
    lai = tmp_path / "lai.npy"
    np.save(green, np.array([[4.9]]))
    np.save(lai, np.array([[0.25]]))
    missing = tmp_path / "unselected.npy"
    for option in (0, 1, 2):
        result = initialize(
            *paths,
            dust_veg=option,
            vegfra_file=green if option == 1 else missing,
            lai_vegmask_file=lai if option == 2 else missing,
        )
        np.testing.assert_array_equal(result[0], values[0])
        if option == 1:
            np.testing.assert_array_equal(result[5], np.array([[4.9]]))
            assert result[6] is None
        elif option == 2:
            assert result[5] is None
            np.testing.assert_array_equal(result[6], np.array([[0.25]]))
        else:
            assert result[5] is None and result[6] is None


def test_a4a5_initialization_requires_selected_dataset(tmp_path):
    """A5: Initialization must fail when the selected vegetation file is absent."""
    paths = [tmp_path / (name + ".npy") for name in ["erod", "clay", "sand", "isltyp", "oro"]]
    values = [
        np.array([[[0.8, 0.0, 0.0]]]),
        np.array([[0.2]]),
        np.array([[0.6]]),
        np.array([[3]]),
        np.array([[1.0]]),
    ]
    for path, value in zip(paths, values):
        np.save(path, value)
    for option, field in [(1, "vegfra_file"), (2, "lai_vegmask_file")]:
        with pytest.raises(ValueError, match=field):
            initialize(*paths, dust_veg=option)


# --- Remaining WRF driver options ---


def test_rest_soil_type_controls_porosity_without_override(cell):
    """New 1: Soil type must determine porosity without an inconsistent caller override."""
    import inspect

    assert "poros" not in inspect.signature(dust_emission_afwa).parameters
    fields = {k: v for k, v in cell.items() if k != "poros"}
    expected = np.array(
        [
            6.58893048854812215e-8,
            6.20856392403230911e-8,
            1.27484148558687142e-7,
            2.95520280850431189e-7,
            6.25150853615506318e-8,
        ]
    )
    np.testing.assert_allclose(dust_emission_afwa(**fields)[0, 0], expected, rtol=1e-12, atol=0)


def test_rest_unknown_vegetation_options_use_ginoux(cell):
    """New 2: Unknown selectors must retain WRF's Ginoux fallback."""
    for erod in [np.array([[[0.8, 0, 0]]]), np.array([[[0, 0.3, 0.5]]])]:
        fields = {**cell, "erod": erod}
        reference = dust_emission_afwa(**fields, dust_veg=0)
        for option in (-1, 3):
            np.testing.assert_array_equal(dust_emission_afwa(**fields, dust_veg=option), reference)


def test_rest_ruc_moisture_uses_volumetric_threshold(cell):
    """A6: RUC must use percent volumetric moisture above its dry limit."""
    fields = {**cell, "volsm": np.array([[0.05]])}
    assert (
        dust_emission_afwa(**fields, dust_smois=1, sf_surface_physics=3).sum()
        < dust_emission_afwa(**fields).sum()
    )


def test_rest_px_moisture_uses_volumetric_threshold(cell):
    """A6: PX must use the RUC/PX volumetric branch with Registry value 7."""
    fields = {**cell, "volsm": np.array([[0.05]])}
    assert (
        dust_emission_afwa(**fields, dust_smois=1, sf_surface_physics=7).sum()
        < dust_emission_afwa(**fields).sum()
    )


def test_rest_noah_moisture_uses_cosby_dry_limit(cell):
    """A6: Noah must use its own volumetric dry limit with Registry value 2."""
    fields = {**cell, "volsm": np.array([[0.055]])}
    assert (
        dust_emission_afwa(**fields, dust_smois=1, sf_surface_physics=2).sum()
        < dust_emission_afwa(**fields).sum()
    )


def test_rest_other_surface_model_retains_gravimetric_moisture(cell):
    """A6: Unsupported surface models must fall back to gravimetric moisture."""
    fields = {**cell, "volsm": np.array([[0.06]])}
    np.testing.assert_array_equal(
        dust_emission_afwa(**fields, dust_smois=1, sf_surface_physics=4),
        dust_emission_afwa(**fields),
    )


def test_rest_other_moisture_selector_retains_gravimetric_moisture(cell):
    """A6: Only dust_smois=1 may select a volumetric branch."""
    np.testing.assert_array_equal(
        dust_emission_afwa(**cell, dust_smois=2, sf_surface_physics=3), dust_emission_afwa(**cell)
    )


def test_rest_dri_present_selects_regional_erodibility(cell):
    """A7: Available DRI data must replace total erodibility but not the Ginoux mask."""
    dri = np.array([[[0.2, 0.2, 0.2]]])
    reference = dust_emission_afwa(**{**cell, "erod": np.array([[[0.6, 0, 0]]])})
    np.testing.assert_allclose(
        dust_emission_afwa(**cell, dust_dsr=1, erod_dri=dri), reference, rtol=1e-14, atol=0
    )
    assert (
        dust_emission_afwa(
            **{**cell, "erod": np.array([[[0, 0.3, 0.5]]])}, dust_dsr=1, erod_dri=dri
        )
        == 0
    ).all()


def test_rest_dri_missing_region_falls_back_to_ginoux(cell):
    """A7: A negative DRI first bin must select the original Ginoux total."""
    np.testing.assert_array_equal(
        dust_emission_afwa(**cell, dust_dsr=1, erod_dri=np.array([[[-1, 0.9, 0.9]]])),
        dust_emission_afwa(**cell),
    )


def test_rest_nga_present_controls_all_texture_physics(cell):
    """A7: NGA textures must reach weighting, beta and both moisture paths."""
    clay, sand = np.array([[0.05]]), np.array([[0.8]])
    for options in [{}, {"dust_smois": 1, "sf_surface_physics": 3}]:
        np.testing.assert_array_equal(
            dust_emission_afwa(**cell, dust_soils=1, clay_nga=clay, sand_nga=sand, **options),
            dust_emission_afwa(**{**cell, "clay": clay, "sand": sand}, **options),
        )


def test_rest_nga_missing_region_falls_back_to_normal_texture(cell):
    """A7: Negative NGA clay must fall back for both clay and sand."""
    np.testing.assert_array_equal(
        dust_emission_afwa(
            **cell, dust_soils=1, clay_nga=np.array([[-1]]), sand_nga=np.array([[0.9]])
        ),
        dust_emission_afwa(**cell),
    )


def test_rest_cmb_cannot_override_wrf_saltation(cell):
    """A8: The source's fixed saltation coefficient must reject caller overrides."""
    with pytest.raises(TypeError, match="cmb"):
        dust_emission_afwa(**cell, cmb=2.61)


def test_rest_selected_regional_options_require_datasets(cell):
    """A7: Selected regional options must not silently replace missing datasets."""
    for options, name in [({"dust_dsr": 1}, "erod_dri"), ({"dust_soils": 1}, "clay_nga")]:
        with pytest.raises(ValueError, match=name):
            dust_emission_afwa(**cell, **options)


def test_rest_initializer_loads_only_selected_regional_data(tmp_path):
    """A7: Static initialization must preserve regional fallback sentinels and ignore unused files."""
    names = ["erod", "clay", "sand", "isltyp", "oro", "erod_dri", "clay_nga", "sand_nga"]
    values = [
        np.array([[[0.8, 0, 0]]]),
        np.array([[0.2]]),
        np.array([[0.6]]),
        np.array([[3]]),
        np.array([[1.0]]),
        np.array([[[-1.0, 0.2, 0.3]]]),
        np.array([[-1.0]]),
        np.array([[0.9]]),
    ]
    paths = [tmp_path / (name + ".npy") for name in names]
    for path, value in zip(paths, values):
        np.save(path, value)
    selected = initialize(
        *paths[:5],
        dust_dsr=1,
        erod_dri_file=paths[5],
        dust_soils=1,
        clay_nga_file=paths[6],
        sand_nga_file=paths[7],
    )
    for actual, expected in zip(selected[7:], values[5:]):
        np.testing.assert_array_equal(actual, expected)
    missing = tmp_path / "unused.npy"
    unused = initialize(
        *paths[:5], dust_veg=-1, erod_dri_file=missing, clay_nga_file=missing, sand_nga_file=missing
    )
    assert all(value is None for value in unused[5:])
    with pytest.raises(ValueError, match="erod_dri_file"):
        initialize(*paths[:5], dust_dsr=1)
    with pytest.raises(ValueError, match="clay_nga_file"):
        initialize(*paths[:5], dust_soils=1)


def test_rest_fractional_soil_type_cannot_bypass_exclusions(cell):
    """New 1: Fractional soil types must not use one class's porosity and another mask."""
    with pytest.raises(ValueError, match="integer"):
        dust_emission_afwa(**{**cell, "isltyp": np.array([[15.5]])})


# --- Review findings: preserve WRF's inactive-cell and comparison semantics ---


def test_review_water_soil_values_cannot_trigger_land_lookup(cell):
    """Water cells must emit zero without validating or looking up their soil class."""
    for soil in [0, -1, 20, 15.5, np.nan]:
        flux = dust_emission_afwa(
            **{**cell, "oro": np.array([[2.0]]), "isltyp": np.array([[soil]])}
        )
        assert (flux == 0).all()


def test_review_initializer_preserves_arbitrary_water_soil(tmp_path):
    """Initialization must not reject a soil value WRF never reads over water."""
    names = ["erod", "clay", "sand", "isltyp", "oro"]
    values = [
        np.array([[[0.8, 0, 0]]]),
        np.array([[0.2]]),
        np.array([[0.6]]),
        np.array([[0.0]]),
        np.array([[2.0]]),
    ]
    paths = [tmp_path / (name + ".npy") for name in names]
    for path, value in zip(paths, values):
        np.save(path, value)
    result = initialize(*paths)
    assert result[3][0, 0] == 0


def test_review_initializer_accepts_wrf_texture_and_erodibility_values(tmp_path):
    """Loading must preserve WRF inputs instead of imposing new physical range checks."""
    names = ["erod", "clay", "sand", "isltyp", "oro"]
    paths = [tmp_path / (name + ".npy") for name in names]
    for erod, clay, sand in [(2.0, 0.6, 0.6), (-0.1, -1.0, 1.2)]:
        values = [
            np.array([[[erod, 0, 0]]]),
            np.array([[clay]]),
            np.array([[sand]]),
            np.array([[3]]),
            np.array([[1.0]]),
        ]
        for path, value in zip(paths, values):
            np.save(path, value)
        actual = initialize(*paths)
        for loaded, expected in zip(actual[:3], values[:3]):
            np.testing.assert_array_equal(loaded, expected)


def test_review_nan_roughness_does_not_satisfy_wrf_suppression(cell):
    """NaN roughness must retain emission because WRF only masks znt above 0.2."""
    np.testing.assert_array_equal(
        dust_emission_afwa(**{**cell, "znt": np.array([[np.nan]])}), dust_emission_afwa(**cell)
    )


def test_review_nan_snow_does_not_satisfy_wrf_suppression(cell):
    """NaN snow must retain emission because WRF only masks snowh above 0.01."""
    np.testing.assert_array_equal(
        dust_emission_afwa(**{**cell, "snowh": np.array([[np.nan]])}), dust_emission_afwa(**cell)
    )


def _driver_test_inputs(cell):
    fields = dict(cell)
    density = fields.pop("airden")
    moisture = fields.pop("volsm")
    return dict(
        dust=np.array([[[[10.0, 20.0, 30.0, 40.0, 50.0]], [[2.0, 4.0, 6.0, 8.0, 10.0]]]]),
        rho_phy=np.stack([density, density * 0.8], axis=1),
        dz8w=np.array([[[50.0], [500.0]]]),
        smois=np.stack([moisture, moisture * 4], axis=1),
        dt=120.0,
        u10=np.array([[10.0]]),
        v10=np.array([[3.0]]),
        emis_dust=np.array([[[1.0, 2.0, 3.0, 4.0, 5.0]]]),
        tot_edust=np.array([[20.0]]),
        dustin=np.array([[[6.0, 7.0, 8.0, 9.0, 10.0]]]),
        **fields,
    )


def test_explicit_default_settings_equal_implicit_defaults(cell):
    np.testing.assert_array_equal(
        dust_emission_afwa(**cell), dust_emission_afwa(**cell, settings=AfwaSettings())
    )


def _scale(value):
    return value * 1.2


def _scale_first(value):
    return (value[0] * 1.2, *value[1:])


# One entry per setting. A setting without an entry fails the test below,
# so a renamed or forgotten case cannot pass silently.
AFWA_SETTING_CHANGES = {
    "alpha": _scale,
    "gamma": _scale,
    "smtune": _scale,
    "ustune": _scale,
    "znt_max": lambda value: 0.005,
    "snow_max": lambda value: -0.001,
    "nonerodible_soil": lambda value: (*value, 3),
    "reff_salt": _scale_first,
    "den_salt": _scale_first,
    "spoint": lambda value: (2, *value[1:]),
    "frac_salt": _scale_first,
    "distr_dust": _scale_first,
    "porosity_table": lambda value: (*value[:2], value[2] * 1.2, *value[3:]),
    "betamax": _scale,
    "cmb": _scale,
    "beta_clay": _scale,
    "beta_offset": _scale,
    "veg_cutoff": lambda value: 6.0,
    "dry_ruc_scale": _scale,
    "dry_ruc_clay": _scale,
    "dry_ruc_offset": _scale,
    "dry_ruc_power": _scale,
    "dry_noah_scale": _scale,
    "dry_noah_clay": _scale,
    "dry_noah_offset": _scale,
    "dry_noah_power": _scale,
    "dry_base_quad": _scale,
    "dry_base_linear": _scale,
    "soil_dens_sand": _scale,
    "soil_dens_clay": _scale,
    "moist_scale": _scale,
    "moist_power": _scale,
    "iw_scale": _scale,
    "iw_cohesion": _scale,
    "iw_drag": _scale,
    "iw_size": _scale,
    "iw_offset": _scale,
    "iw_size_power": _scale,
    "iw_drag_power": _scale,
    "iw_diam_power": _scale,
    # Scaling it cancels in the normalized area fraction, so 0 (0/0 = NaN) is the probe.
    "surface_geometry": lambda value: 0.0,
    "loft_height": _scale,
    "extinction": _scale_first,
    "visibility_scale": _scale,
}


def _afwa_outputs(settings):
    wet = {
        **_afwa_cell(),
        "ustar": np.array([[1.5]]),
        "volsm": np.array([[0.12]]),
        "clay": np.array([[0.1]]),
    }
    cases = [
        wet,
        {**wet, "clay": np.array([[0.3]])},
        {**wet, "dust_smois": 1, "sf_surface_physics": 3},
        {**wet, "dust_smois": 1, "sf_surface_physics": 2},
        {**wet, "dust_veg": 1, "vegfra": np.array([[5.5]])},
    ]
    result = [dust_emission_afwa(**sample, settings=settings).ravel() for sample in cases]
    driver = afwa_dust_driver(**_driver_test_inputs(wet), settings=settings)
    result.extend(driver[name].ravel() for name in ("afwa_dustloft", "vis_dust"))
    calm = {**wet, "ustar": np.array([[0.0]])}
    result.append(
        afwa_dust_driver(**_driver_test_inputs(calm), settings=settings)["afwa_dustloft"].ravel()
    )
    return np.concatenate(result)


@pytest.mark.parametrize("name", tunable_names(AfwaSettings))
def test_each_setting_reaches_an_output_so_tuning_it_is_never_silently_ignored(name):
    base = AfwaSettings()
    changed = AFWA_SETTING_CHANGES[name](getattr(base, name))
    assert not np.array_equal(
        _afwa_outputs(replace(base, **{name: changed})), _afwa_outputs(base), equal_nan=True
    ), name


def test_every_setting_case_names_a_real_setting():
    assert set(AFWA_SETTING_CHANGES) <= set(tunable_names(AfwaSettings))


def test_loft_height_enters_both_log_terms_of_the_lofting_potential(cell):
    """The reference height must set the log term and psi, or a tuned height is half applied."""
    height = 20.0
    inputs = _driver_test_inputs(cell)
    inputs.update(
        settings=AfwaSettings(ustune=0.7, loft_height=height), smois=np.array([[[0.2], [0.8]]])
    )
    result = afwa_dust_driver(**inputs)
    clay = cell["clay"]
    moisture = (
        100
        * inputs["smois"][:, 0, :]
        / ((1 - POROSITY_TABLE[2]) * (2.65 * (1 - clay) + 2.5 * clay))
    )
    drylimit = 14 * clay**2 + 17 * clay
    threshold = threshold_friction_velocity(REFF_SALT, DEN_SALT, cell["airden"])[..., 6]
    threshold *= np.sqrt(1 + 1.21 * (moisture - drylimit) ** 0.68)
    wind = np.sqrt(inputs["u10"] ** 2 + inputs["v10"] ** 2)
    psi = 0.4 * wind / (cell["ustar"] * 0.7) - np.log(height / cell["znt"])
    expected = 0.7 * wind - threshold * (np.log(height / cell["znt"]) + psi) / 0.4
    np.testing.assert_allclose(result["afwa_dustloft"], expected, rtol=1e-14)


def test_default_emission_and_driver_numbers_are_pinned_so_a_default_drift_fails(cell):
    """Defaults are the Fortran values; a 1e-6 drift in any of them must show here."""
    wet = {
        **cell,
        "ustar": np.array([[1.5]]),
        "volsm": np.array([[0.12]]),
        "clay": np.array([[0.1]]),
    }
    pinned = [
        (
            cell,
            [
                6.588930488548123e-08,
                6.20856392403231e-08,
                1.2748414855868717e-07,
                2.9552028085043124e-07,
                6.251508536155063e-08,
            ],
        ),
        (
            wet,
            [
                9.308717984297907e-08,
                8.771343203081455e-08,
                1.8010722505932077e-07,
                4.1750553566446017e-07,
                8.832014549347828e-08,
            ],
        ),
        (
            {**wet, "dust_smois": 1, "sf_surface_physics": 3},
            [
                7.486789869226266e-08,
                7.054591571375215e-08,
                1.448561391829812e-07,
                3.3579019366911477e-07,
                7.103388153390656e-08,
            ],
        ),
        (
            {**wet, "dust_smois": 1, "sf_surface_physics": 2},
            [
                7.531725783413496e-08,
                7.096933419752754e-08,
                1.4572556962693897e-07,
                3.378056154441602e-07,
                7.146022880166063e-08,
            ],
        ),
    ]
    for sample, expected in pinned:
        flux = dust_emission_afwa(**sample).ravel()
        assert flux.min() > 0
        np.testing.assert_allclose(flux, expected, rtol=1e-9, atol=0)
    driver = afwa_dust_driver(**_driver_test_inputs(wet))
    np.testing.assert_allclose(
        driver["afwa_dustloft"].ravel(), [7.165410841688542], rtol=1e-9, atol=0
    )
    np.testing.assert_allclose(
        driver["vis_dust"].ravel(), [3755.9521511204407, 327339.2540012035], rtol=1e-9, atol=0
    )


def test_driver_dust_updates_only_lowest_level_in_ug_per_kg(cell):
    """Emitted mass must increase the lowest tracers without changing higher levels."""
    inputs = _driver_test_inputs(cell)
    output = afwa_dust_driver(**inputs)
    flux = dust_emission_afwa(**cell)
    expected = (
        inputs["dust"][:, 0, :, :] * 1.0e-9
        + flux * inputs["dt"] / inputs["dz8w"][:, 0, :, None] / cell["airden"][..., None]
    ) * 1.0e9
    np.testing.assert_array_equal(output["dust"][:, 0, :, :], expected)
    np.testing.assert_array_equal(output["dust"][:, 1, :, :], inputs["dust"][:, 1, :, :])
    np.testing.assert_array_equal(inputs["dust"][0, 0, 0], [10.0, 20.0, 30.0, 40.0, 50.0])


def test_driver_loft_uses_bin_seven_moisture_threshold(cell):
    """Lofting potential must use the wet small-sand threshold, even below emission cutoff."""
    inputs = _driver_test_inputs(cell)
    inputs.update(settings=AfwaSettings(ustune=0.7), smois=np.array([[[0.2], [0.8]]]))
    result = afwa_dust_driver(**inputs)
    clay = cell["clay"]
    moisture = (
        100
        * inputs["smois"][:, 0, :]
        / ((1 - POROSITY_TABLE[2]) * (2.65 * (1 - clay) + 2.5 * clay))
    )
    drylimit = 14 * clay**2 + 17 * clay
    threshold = threshold_friction_velocity(REFF_SALT, DEN_SALT, cell["airden"])[..., 6]
    threshold *= np.sqrt(1 + 1.21 * (moisture - drylimit) ** 0.68)
    wind = np.sqrt(inputs["u10"] ** 2 + inputs["v10"] ** 2)
    psi = 0.4 * wind / (cell["ustar"] * 0.7) - np.log(10 / cell["znt"])
    expected = 0.7 * wind - threshold * (np.log(10 / cell["znt"]) + psi) / 0.4
    np.testing.assert_allclose(result["afwa_dustloft"], expected, rtol=1e-14)


def test_driver_emis_dust_accumulates_mass_not_flux(cell):
    """Each emitted-bin total must retain its initial mass and add flux times dt."""
    inputs = _driver_test_inputs(cell)
    result = afwa_dust_driver(**inputs)
    np.testing.assert_array_equal(
        result["emis_dust"], inputs["emis_dust"] + dust_emission_afwa(**cell) * 120.0
    )


def test_driver_tot_edust_retains_running_total(cell):
    """Total emitted dust must add all five bin masses to the existing accumulator."""
    inputs = _driver_test_inputs(cell)
    result = afwa_dust_driver(**inputs)
    np.testing.assert_array_equal(
        result["tot_edust"], inputs["tot_edust"] + (dust_emission_afwa(**cell) * 120.0).sum(axis=-1)
    )


def test_driver_tot_dust_uses_all_levels_and_density(cell):
    """Concentration diagnostics must use updated chemistry and density at every level."""
    inputs = _driver_test_inputs(cell)
    result = afwa_dust_driver(**inputs)
    np.testing.assert_array_equal(
        result["tot_dust"], result["dust"].sum(axis=-1) * inputs["rho_phy"]
    )


def test_driver_visibility_uses_extinction_and_actual_wrf_cap(cell):
    """Visibility must follow the five-bin extinction coefficients and 999999 m cap."""
    inputs = _driver_test_inputs(cell)
    result = afwa_dust_driver(**inputs)
    coeff = np.array([1.470e-6, 7.877e-7, 4.623e-7, 2.429e-7, 1.387e-7])
    expected = np.minimum(
        3.912 / ((result["dust"] * coeff).sum(axis=-1) * inputs["rho_phy"]), 999999.0
    )
    np.testing.assert_allclose(result["vis_dust"], expected, rtol=1e-14)
    inputs["dust"] = np.zeros_like(inputs["dust"])
    inputs["ustar"] = np.zeros_like(inputs["ustar"])
    np.testing.assert_array_equal(afwa_dust_driver(**inputs)["vis_dust"], 999999.0)


@pytest.mark.parametrize("chem_opt", [2, 11])
def test_driver_dustin_routes_updated_tracers_without_changing_chemistry(cell, chem_opt):
    """Chemistry options 2 and 11 update dustin while concentration uses unchanged chem."""
    inputs = _driver_test_inputs(cell)
    result = afwa_dust_driver(**inputs, chem_opt=chem_opt)
    expected = (
        inputs["dust"][:, 0, :, :] * 1.0e-9
        + dust_emission_afwa(**cell) * 120.0 / 50.0 / cell["airden"][..., None]
    ) * 1.0e9
    np.testing.assert_array_equal(result["dustin"], expected)
    np.testing.assert_array_equal(result["dust"], inputs["dust"])
    np.testing.assert_array_equal(
        result["tot_dust"], inputs["dust"].sum(axis=-1) * inputs["rho_phy"]
    )


@pytest.mark.parametrize("density", [0.0, -1.0])
def test_driver_sea_skips_invalid_texture_density_without_warnings(cell, density):
    """Sea emission work must be skipped while retaining WRF's sea diagnostics and state."""
    import warnings

    sea = dict(
        cell,
        oro=np.array([[2.0]]),
        isltyp=np.array([[0]]),
        clay=np.array([[np.nan]]),
        sand=np.array([[np.inf]]),
        airden=np.array([[density]]),
    )
    inputs = _driver_test_inputs(sea)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        flux = dust_emission_afwa(**sea)
        result = afwa_dust_driver(**inputs)
    np.testing.assert_array_equal(flux, 0.0)
    for name in ("dust", "emis_dust", "tot_edust", "dustin"):
        np.testing.assert_array_equal(result[name], inputs[name])
    np.testing.assert_array_equal(result["afwa_dustloft"], -99.0)
    np.testing.assert_array_equal(
        result["tot_dust"], inputs["dust"].sum(axis=-1) * inputs["rho_phy"]
    )
    np.testing.assert_array_equal(result["vis_dust"], 999999.0)
    assert all(np.isfinite(value).all() for value in result.values())


def test_driver_mixed_grid_keeps_land_tuning_separate_from_invalid_sea(cell):
    """Packing land cells must preserve both horizontal axes and each cell's dt/tuning."""
    import warnings

    inputs = _driver_test_inputs(cell)
    for name, value in list(inputs.items()):
        if np.ndim(value) == 0:
            continue
        shape = (
            (2, 2, 2, 5)
            if name == "dust"
            else (2, 2, 2)
            if name in ("rho_phy", "dz8w", "smois")
            else (2, 2) + value.shape[2:]
        )
        inputs[name] = np.broadcast_to(value, shape).copy()
    sea = np.array([[False, True], [True, False]])
    inputs["oro"][sea] = 2.0
    inputs["isltyp"][sea] = 0
    inputs["clay"][sea] = np.nan
    inputs["sand"][sea] = np.inf
    inputs["rho_phy"][:, 0, :][sea] = [0.0, -1.0]
    inputs["settings"] = AfwaSettings(ustune=np.array([[1.0, np.nan], [np.inf, 0.7]]))
    inputs["dt"] = np.array([[120.0, np.nan], [np.inf, 60.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = afwa_dust_driver(**inputs)
    for i, j in ((0, 0), (1, 1)):
        single = _driver_test_inputs(cell)
        single.update(
            dt=inputs["dt"][i, j], settings=AfwaSettings(ustune=inputs["settings"].ustune[i, j])
        )
        expected = afwa_dust_driver(**single)
        for name in result:
            actual = (
                result[name][i, :, j, ...]
                if name in ("dust", "tot_dust", "vis_dust")
                else result[name][i, j, ...]
            )
            target = (
                expected[name][0, :, 0, ...]
                if name in ("dust", "tot_dust", "vis_dust")
                else expected[name][0, 0, ...]
            )
            np.testing.assert_array_equal(actual, target)
    np.testing.assert_array_equal(result["flux"][sea], 0.0)
    np.testing.assert_array_equal(result["afwa_dustloft"][sea], -99.0)


def test_driver_uses_lowest_dz8w_level(cell):
    inputs = _driver_test_inputs(cell)
    base = afwa_dust_driver(**inputs)
    changed = dict(inputs, dz8w=inputs["dz8w"].copy())
    changed["dz8w"][:, 1, :] = 5.0
    same = afwa_dust_driver(**changed)
    np.testing.assert_array_equal(same["dust"], base["dust"])
    changed["dz8w"][:, 0, :] = 100.0
    assert not np.array_equal(afwa_dust_driver(**changed)["dust"], base["dust"])


def test_driver_uses_first_soil_layer(cell):
    inputs = _driver_test_inputs(cell)
    base = afwa_dust_driver(**inputs)
    changed = dict(inputs, smois=inputs["smois"].copy())
    changed["smois"][:, 1, :] = 0.9
    np.testing.assert_array_equal(afwa_dust_driver(**changed)["flux"], base["flux"])
    changed["smois"][:, 0, :] = 0.2
    assert not np.array_equal(afwa_dust_driver(**changed)["flux"], base["flux"])
