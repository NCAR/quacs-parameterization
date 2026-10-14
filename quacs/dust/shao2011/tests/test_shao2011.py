"""Tests for quacs.dust.shao2011."""

import warnings
from dataclasses import replace

import numpy as np
import pytest

import quacs.dust.shao2011.shao2011 as shao
from quacs.dust.shao2011.shao2011 import (
    C_WHITE,
    CC2ISL,
    CYS_DEFAULT,
    DBIN,
    F_SHAO,
    GAMMA_SHAO,
    GRAV,
    IMAX,
    N_BINS,
    N_SOIL,
    NOAH_DRYSMC,
    PHL,
    RHOP_DEFAULT,
    RHOS,
    RUC_DRYSMC,
    SOIL_CLASS_NAMES,
    SOIL_PSD_F,
    SOIL_PSD_M,
    WR_THR,
    Shao2011Settings,
    ShaoDustWarning,
    _flux_and_diagnostics,
    dust_emission_shao2011,
    initialize,
    moisture_correction_fecan,
    psd_create,
    roughness_correction_raupach,
    saltation_flux_white,
    threshold_velocity_shao_lu,
    tunable_names,
    uoc_dust_driver,
)

# ===========================================================================
# Pytest test suite
# Run with: pytest dust_shao2011.py -v
# ===========================================================================


@pytest.fixture(scope="module")
def psd():
    return initialize()


@pytest.fixture()
def single_cell():
    """Single active sandy cell, dry soil, strong wind, sparse cover."""
    soil = np.zeros((1, 1, N_SOIL))
    soil[0, 0, 0] = 1.0  # 100% sand
    return dict(
        ustar=np.array([[0.8]]),
        gwet=np.array([[0.02]]),
        cf=np.array([[5.0]]),
        erod=np.array([[[0.5, 0.0, 0.0]]]),
        soil_frac=soil,
        rhoa=np.array([[1.225]]),
    )


def test_explicit_default_settings_equal_implicit_defaults(single_cell, psd):
    np.testing.assert_array_equal(
        dust_emission_shao2011(**single_cell, psd=psd),
        dust_emission_shao2011(**single_cell, psd=psd, settings=Shao2011Settings()),
    )


def _scale(value):
    return value * 1.2


def _scale_first(value):
    return (value[0] * 1.2, *value[1:])


def _scale_first_texture(value):
    changed = list(value)
    first = list(changed[0])
    first[0] = (first[0][0] * 1.2, *first[0][1:])
    changed[0] = tuple(first)
    return tuple(changed)


# One entry per setting. A setting without an entry fails the test below.
SHAO_SETTING_CHANGES = {
    "c_y": _scale,
    "rhos": _scale,
    "phl": _scale,
    "rhop": _scale,
    "f_shao": _scale,
    "gamma_shao": _scale,
    "c_white": _scale,
    "c_lambda": _scale,
    "raup_sig": _scale,
    "raup_m": _scale,
    "raup_beta": _scale,
    "cover_cap": lambda value: 0.9,
    "imax": lambda value: value + 10,
    "jmax": lambda value: value - 1,
    "dref": _scale,
    "phi_up": _scale,
    "phi_lo": _scale,
    "dcut": _scale,
    "dbin": lambda value: (3.0, *value[1:]),
    "wr_thr": _scale_first,
    "noah_drysmc": _scale_first,
    "ruc_drysmc": _scale_first,
    "fecan_a": _scale_first,
    "fecan_b": _scale_first,
    "soil_psd_m": _scale_first_texture,
    "soil_psd_f": _scale_first_texture,
    "ghl_power": _scale,
    "sigma_scale": _scale,
    "sigma_slope": _scale,
    "impact_angle": _scale,
    "impact_wind": _scale,
    "impact_cubic": _scale,
}


def _shao_case(single_cell):
    return {**single_cell, "ustar": np.array([[1.5]]), "gwet": np.array([[0.08]])}


def _shao_outputs(settings, single_cell):
    case = _shao_case(single_cell)
    psd_value = initialize(settings=settings)
    result = [
        psd_value["dm"].ravel(),
        psd_value["dpsd_m"].ravel(),
        psd_value["dpsd_f"].ravel(),
        psd_value["bin_hi"].ravel(),
        np.array([psd_value["idst"]]),
        roughness_correction_raupach(np.array([[0.97]]), settings).ravel(),
    ]
    for sf in (2, 3, 0):
        for imod in (1, 2, 3):
            result.append(
                dust_emission_shao2011(
                    **case, psd=psd_value, sf_surface_physics=sf, imod=imod, settings=settings
                ).ravel()
            )
    return np.concatenate(result)


@pytest.mark.parametrize("name", tunable_names(Shao2011Settings))
def test_each_setting_reaches_an_output_so_tuning_it_is_never_silently_ignored(name, single_cell):
    base = Shao2011Settings()
    changed = SHAO_SETTING_CHANGES[name](getattr(base, name))
    assert not np.array_equal(
        _shao_outputs(replace(base, **{name: changed}), single_cell),
        _shao_outputs(base, single_cell),
        equal_nan=True,
    ), name


def test_every_setting_case_names_a_real_setting():
    assert set(SHAO_SETTING_CHANGES) <= set(tunable_names(Shao2011Settings))


def test_particle_density_enters_the_imod_1_impact_mass_not_only_the_threshold(single_cell):
    """rhop feeds the dry threshold and the imod 1 impact mass; the mass use must read the setting."""
    case = _shao_case(single_cell)
    base, k = Shao2011Settings(), 1.5
    # Scaling f_shao by 1/k and gamma_shao by k keeps the dry threshold fixed,
    # so any change in the flux comes from the impact mass alone.
    tuned = replace(
        base, rhop=base.rhop * k, f_shao=base.f_shao / k, gamma_shao=base.gamma_shao * k
    )
    same = dust_emission_shao2011(**case, psd=initialize(), imod=3)
    other = dust_emission_shao2011(**case, psd=initialize(settings=tuned), imod=3, settings=tuned)
    np.testing.assert_allclose(other, same, rtol=1e-12)
    flux = dust_emission_shao2011(**case, psd=initialize(), imod=1)
    tuned_flux = dust_emission_shao2011(
        **case, psd=initialize(settings=tuned), imod=1, settings=tuned
    )
    assert (np.abs(tuned_flux - flux) / flux > 0.05).all()


def test_default_numbers_are_pinned_so_a_default_drift_fails(single_cell):
    """Defaults are the WRF values; a 1e-6 drift in any of them must show here."""
    case = _shao_case(single_cell)
    psd_value = initialize()
    pinned = {
        (2, 1): [
            4.234159230439125e-42,
            3.352466798610835e-29,
            1.8703373691383983e-21,
            6.486076524558259e-13,
            8.690375871877219e-09,
        ],
        (2, 2): [
            7.126368070190844e-42,
            5.642421800825247e-29,
            3.1479006297385527e-21,
            1.0916492774560268e-12,
            1.4626473347078134e-08,
        ],
        (2, 3): [
            8.440700452201393e-31,
            8.3896846270819615e-22,
            2.0257825992349548e-16,
            1.6121605750476955e-10,
            1.0828906547166185e-07,
        ],
        (3, 1): [
            4.185461638344465e-42,
            3.313909660868826e-29,
            1.8488263863612456e-21,
            6.411479351495657e-13,
            8.590426777777124e-09,
        ],
        (3, 2): [
            7.044406824401188e-42,
            5.577527605702098e-29,
            3.11169623295548e-21,
            1.0790940833004048e-12,
            1.445825245738676e-08,
        ],
        (3, 3): [
            8.34359794178124e-31,
            8.2931690068989745e-22,
            2.00247783003184e-16,
            1.5936141475415867e-10,
            1.0704329917917265e-07,
        ],
    }
    for (sf, imod), expected in pinned.items():
        flux = dust_emission_shao2011(
            **case, psd=psd_value, sf_surface_physics=sf, imod=imod
        ).ravel()
        assert flux.min() > 0
        np.testing.assert_allclose(flux, expected, rtol=1e-9, atol=0, err_msg=str((sf, imod)))


# --- source-pinned constants ---


def test_constants_pinned_to_source():
    assert CYS_DEFAULT == 1.0e-5
    assert RHOS == 1000.0
    assert PHL == 30000.0
    assert RHOP_DEFAULT == 2650.0
    assert F_SHAO == 0.0123
    assert GAMMA_SHAO == 1.65e-4
    assert C_WHITE == 2.3


def test_bin_cut_diameters():
    assert np.allclose(DBIN, [2.0, 3.6, 6.0, 12.0, 20.0])


def test_soil_tables_complete():
    assert len(SOIL_PSD_M) == N_SOIL
    assert len(SOIL_PSD_F) == N_SOIL
    for name in SOIL_CLASS_NAMES:
        assert name in SOIL_PSD_M and name in SOIL_PSD_F


# --- PSD construction ---


def test_psd_normalized():
    """psd_create must return a normalized cumulative PSD (P(dmax) = 1)."""
    _, _, _, dpsd, ppsd = psd_create(SOIL_PSD_M["sand"])
    assert np.isclose(ppsd[-1], 1.0, rtol=1e-10)
    assert np.isclose(dpsd.sum(), 1.0, rtol=1e-10)
    assert (dpsd >= -1e-12).all()


def test_psd_uses_source_pi_not_numpy_pi(monkeypatch):
    """S7: psd_create's Fortran 3.14159 is independent of NumPy's pi."""
    original = psd_create(SOIL_PSD_M["sand"])
    monkeypatch.setattr(np, "pi", 3.0)
    changed = psd_create(SOIL_PSD_M["sand"])
    assert all(np.array_equal(a, b) for a, b in zip(original, changed))


def test_size_grid_monotonic():
    d, dm, _, _, _ = psd_create(SOIL_PSD_M["clay"])
    assert (np.diff(d) > 0).all()
    assert (np.diff(dm) > 0).all()
    assert dm[0] < dm[-1]


def test_initialize_shapes(psd):
    assert psd["dm"].shape == (IMAX,)
    assert psd["d"].shape == (IMAX + 1,)
    assert psd["dpsd_m"].shape == (N_SOIL, IMAX)
    assert psd["dpsd_f"].shape == (N_SOIL, IMAX)
    assert 0 < psd["idst"] <= IMAX
    assert psd["bin_hi"].shape == (N_BINS,)


def test_ibin_search_stops_when_dust_class_is_missing(monkeypatch):
    """Q:210-216: missing upper bin boundaries are an error, not a clamp."""
    monkeypatch.setattr(shao, "DBIN", np.array([2.0, 3.6, 6.0, 12.0, 1.0e9]))
    with pytest.raises(ValueError, match="wrong dust classes"):
        initialize()


# --- threshold velocity (Shao & Lu 2000) ---


def test_threshold_has_interior_minimum():
    """
    Shao-Lu u*t0(d) is U-shaped: cohesion dominates fine particles, gravity
    dominates coarse ones, with a minimum at an intermediate (sand) size.
    """
    d = np.logspace(0, 3, 400)  # 1 .. 1000 micron
    sigma = RHOP_DEFAULT / 1.225
    u = threshold_velocity_shao_lu(d, sigma, 1.225)
    imin = u.argmin()
    assert 0 < imin < len(d) - 1, "threshold minimum should be interior"
    d_min = d[imin]
    assert 30.0 < d_min < 150.0, f"minimum at {d_min:.1f} um, expected sand range"


def test_threshold_increases_with_moisture():
    """Fecan correction can only raise the threshold (H >= 1)."""
    h_dry = moisture_correction_fecan(np.array([0.01]), 0.05, 0)
    h_wet = moisture_correction_fecan(np.array([0.20]), 0.05, 0)
    assert np.isclose(h_dry[0], 1.0)
    assert h_wet[0] > 1.0


# --- saltation flux (White 1979) ---


def test_saltation_zero_below_threshold():
    q = saltation_flux_white(np.array([0.6]), np.array([0.3]), np.array([1.225]))
    assert q[0] == 0.0


def test_saltation_increases_with_ustar():
    u_t = np.array([0.3])
    q_low = saltation_flux_white(u_t, np.array([0.4]), np.array([1.225]))
    q_high = saltation_flux_white(u_t, np.array([0.8]), np.array([1.225]))
    assert q_high[0] > q_low[0] > 0.0


# --- core emission contract ---


def test_output_shape(single_cell, psd):
    F = dust_emission_shao2011(psd=psd, **single_cell)
    assert F.shape == (1, 1, N_BINS)


def test_all_fluxes_nonnegative(single_cell, psd):
    F = dust_emission_shao2011(psd=psd, **single_cell)
    assert (F >= 0.0).all()


def test_active_cell_emits(single_cell, psd):
    F = dust_emission_shao2011(psd=psd, **single_cell)
    assert F.sum() > 0.0


def test_zero_emission_where_erodibility_zero(single_cell, psd):
    """Erodibility is a location gate: erod=0 -> no emission, regardless of wind."""
    c = dict(single_cell)
    c["erod"] = np.zeros((1, 1, 3))
    F = dust_emission_shao2011(psd=psd, **c)
    assert (F == 0.0).all()


def test_erodibility_magnitude_does_not_scale_flux(single_cell, psd):
    """
    Unlike Ginoux, the erodibility VALUE must not scale the flux. Two active
    cells with erod=0.1 and erod=0.9 (all else equal) must emit identically.
    """
    c_lo = dict(single_cell)
    c_lo["erod"] = np.array([[[0.1, 0.0, 0.0]]])
    c_hi = dict(single_cell)
    c_hi["erod"] = np.array([[[0.9, 0.0, 0.0]]])
    F_lo = dust_emission_shao2011(psd=psd, **c_lo)
    F_hi = dust_emission_shao2011(psd=psd, **c_hi)
    assert np.allclose(F_lo, F_hi)


def test_zero_emission_below_threshold(single_cell, psd):
    """Very low u* (below all sand thresholds) gives zero emission."""
    c = dict(single_cell)
    c["ustar"] = np.array([[0.05]])
    F = dust_emission_shao2011(psd=psd, **c)
    assert (F == 0.0).all()


def test_wet_soil_suppresses_emission(single_cell, psd):
    """Higher soil moisture raises the threshold and reduces emission."""
    c_dry = dict(single_cell)
    c_dry["gwet"] = np.array([[0.02]])
    c_wet = dict(single_cell)
    c_wet["gwet"] = np.array([[0.30]])
    F_dry = dust_emission_shao2011(psd=psd, **c_dry)
    F_wet = dust_emission_shao2011(psd=psd, **c_wet)
    assert F_wet.sum() < F_dry.sum()


def test_emission_increases_with_ustar(single_cell, psd):
    c_lo = dict(single_cell)
    c_lo["ustar"] = np.array([[0.5]])
    c_hi = dict(single_cell)
    c_hi["ustar"] = np.array([[1.0]])
    F_lo = dust_emission_shao2011(psd=psd, **c_lo)
    F_hi = dust_emission_shao2011(psd=psd, **c_hi)
    assert F_hi.sum() > F_lo.sum()


def test_soil_fraction_weighting(single_cell, psd):
    """
    A 50/50 sand+clay cell must emit a soil-area-weighted average of the two
    pure-soil emissions (linear in the soil fractions for fixed forcing).
    """
    base = dict(single_cell)
    pure_sand = np.zeros((1, 1, N_SOIL))
    pure_sand[0, 0, 0] = 1.0
    pure_clay = np.zeros((1, 1, N_SOIL))
    pure_clay[0, 0, 11] = 1.0
    mix = np.zeros((1, 1, N_SOIL))
    mix[0, 0, 0] = 0.5
    mix[0, 0, 11] = 0.5

    F_sand = dust_emission_shao2011(psd=psd, **{**base, "soil_frac": pure_sand})
    F_clay = dust_emission_shao2011(psd=psd, **{**base, "soil_frac": pure_clay})
    F_mix = dust_emission_shao2011(psd=psd, **{**base, "soil_frac": mix})
    assert np.allclose(F_mix, 0.5 * F_sand + 0.5 * F_clay, rtol=1e-10)


def test_full_cover_suppresses_emission(single_cell, psd):
    """Dense vegetation (cf -> high) should strongly suppress emission."""
    c = dict(single_cell)
    c["cf"] = np.array([[90.0]])
    F = dust_emission_shao2011(psd=psd, **c)
    F_open = dust_emission_shao2011(psd=psd, **single_cell)
    assert F.sum() < F_open.sum()


def test_grid_vectorization(psd):
    """A 2x2 grid with mixed activity must emit only in the active cell."""
    soil = np.zeros((2, 2, N_SOIL))
    soil[..., 0] = 1.0
    F = dust_emission_shao2011(
        ustar=np.array([[0.8, 0.8], [0.05, 0.8]]),
        gwet=np.full((2, 2), 0.02),
        cf=np.full((2, 2), 5.0),
        erod=np.array([[[0.5, 0.0, 0.0], [0.0, 0.0, 0.0]], [[0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]]),
        soil_frac=soil,
        rhoa=np.full((2, 2), 1.225),
        psd=psd,
    )
    assert F[0, 0].sum() > 0.0  # active, windy
    assert F[0, 1].sum() == 0.0  # erodibility gate off
    assert F[1, 0].sum() == 0.0  # below threshold
    assert F[1, 1].sum() > 0.0  # active, windy


def test_wrf_gravity_and_dry_tables_are_selected(single_cell, psd):
    """S1/S2: gravity and host DRYSMC change the moisture threshold."""
    assert GRAV == 9.81
    assert np.array_equal(
        WR_THR, [0.001, 0.003, 0.037, 0.061, 0.072, 0.049, 0.084, 0.110, 0.095, 0.126, 0.141, 0.156]
    )
    assert np.array_equal(
        NOAH_DRYSMC,
        [0.010, 0.028, 0.047, 0.084, 0.061, 0.066, 0.069, 0.120, 0.103, 0.100, 0.126, 0.138],
    )
    assert np.array_equal(
        RUC_DRYSMC,
        [0.002, 0.035, 0.041, 0.034, 0.034, 0.050, 0.068, 0.060, 0.050, 0.070, 0.070, 0.068],
    )
    c = {**single_cell, "gwet": np.array([[0.02]]), "ustar": np.array([[1.0]])}
    noah = dust_emission_shao2011(psd=psd, **c)
    fallback = dust_emission_shao2011(psd=psd, sf_surface_physics=0, **c)
    assert not np.array_equal(noah, fallback)
    assert np.array_equal(noah, dust_emission_shao2011(psd=psd, sf_surface_physics=4, **c))
    assert not np.array_equal(noah, dust_emission_shao2011(psd=psd, sf_surface_physics=3, **c))


def test_moisture_above_one_clamps_and_negative_active_errors(single_cell, psd):
    """S3: WRF caps active wet cells and h_c stops on negative moisture."""
    wet = {**single_cell, "gwet": np.array([[1.5]])}
    with pytest.warns(ShaoDustWarning, match="moisture"):
        capped = dust_emission_shao2011(psd=psd, **wet)
    assert np.array_equal(
        capped, dust_emission_shao2011(psd=psd, **{**wet, "gwet": np.array([[1.0]])})
    )
    with pytest.raises(ValueError, match="h_c"):
        dust_emission_shao2011(psd=psd, **{**wet, "gwet": np.array([[-0.1]])})
    dust_emission_shao2011(
        psd=psd, **{**wet, "gwet": np.array([[-0.1]]), "erod": np.zeros((1, 1, 3))}
    )


def test_negative_soil_fraction_enters_wrf_weight(single_cell, psd):
    """S4: negative nonzero fractions remain in WRF's normalized sum."""
    soil = np.zeros((1, 1, 12))
    soil[0, 0, :2] = [1.25, -0.25]
    with pytest.warns(ShaoDustWarning, match="fraction"):
        mixed = dust_emission_shao2011(psd=psd, **{**single_cell, "soil_frac": soil})
    sand = dust_emission_shao2011(psd=psd, **single_cell)
    other = soil.copy()
    other[0, 0, :2] = [0.0, 1.0]
    loamysand = dust_emission_shao2011(psd=psd, **{**single_cell, "soil_frac": other})
    assert np.allclose(mixed, 1.25 * sand - 0.25 * loamysand, rtol=1e-12, atol=0)


def test_mixed_sign_soil_sum_uses_wrf_class_order(single_cell, psd):
    """U:209-212: ordered addition preserves a near-zero mixed-sign total."""
    soil = np.zeros((1, 1, 12))
    soil[0, 0, :4] = [0.58, 0.4, -0.33, -0.65]
    assert soil.sum(axis=-1)[0, 0] == 0.0
    with pytest.warns(ShaoDustWarning):
        rate = dust_emission_shao2011(psd=psd, **{**single_cell, "soil_frac": soil})
    assert np.isfinite(rate).all() and np.any(rate != 0.0)


def test_nonfinite_threshold_takes_wrf_emission_branch(single_cell, psd):
    """Q:188: a NaN threshold does not satisfy Fortran's <= no-erosion test."""
    for imod in (1, 2, 3):
        with pytest.warns(ShaoDustWarning):
            rate = dust_emission_shao2011(
                psd=psd, imod=imod, **{**single_cell, "cf": np.array([[np.nan]])}
            )
        assert np.isnan(rate).any()


def test_negative_moisture_with_no_soil_class_does_not_call_h_c(single_cell, psd):
    """Q:116,168: h_c is skipped when every class fraction is zero."""
    soil = np.zeros((1, 1, 12))
    with pytest.warns(ShaoDustWarning, match="fractions are zero"):
        rate = dust_emission_shao2011(
            psd=psd, **{**single_cell, "gwet": np.array([[-0.1]]), "soil_frac": soil}
        )
    assert np.array_equal(rate, np.zeros((1, 1, 5)))


def test_cover_percent_and_three_erod_parts(single_cell, psd):
    """S10/S11: cover uses percent while erod uses only the three-part sum."""
    baseline = dust_emission_shao2011(psd=psd, **single_cell)
    changed = {**single_cell, "erod": np.array([[[0.0, -0.1, 0.6]]])}
    with pytest.warns(ShaoDustWarning, match="erod"):
        assert np.array_equal(baseline, dust_emission_shao2011(psd=psd, **changed))
    with pytest.warns(ShaoDustWarning, match="cover"):
        over = dust_emission_shao2011(psd=psd, **{**single_cell, "cf": np.array([[110.0]])})
    at_cap = dust_emission_shao2011(psd=psd, **{**single_cell, "cf": np.array([[95.0]])})
    assert np.array_equal(over, at_cap)
    with pytest.raises(ValueError, match="erod"):
        dust_emission_shao2011(psd=psd, **{**single_cell, "erod": np.array([[0.5]])})


def test_all_three_imod_branches_and_validation(single_cell, psd):
    """Q:223-339: all three schemes must produce finite distinct bin rates."""
    rates = [dust_emission_shao2011(psd=psd, imod=i, **single_cell) for i in (1, 2, 3)]
    assert all(np.isfinite(rate).all() and rate.sum() > 0 for rate in rates)
    assert not np.array_equal(rates[0], rates[1])
    assert not np.array_equal(rates[1], rates[2])
    with pytest.raises(ValueError, match="imod"):
        dust_emission_shao2011(psd=psd, imod=4, **single_cell)


def test_driver_updates_tracers_totals_moisture_and_diagnostics(single_cell, psd):
    """U:183-266: driver applies mass, conversion, clamp and diagnostics."""
    soil = np.moveaxis(single_cell["soil_frac"], -1, 1)
    start = np.full((1, 1, 1, 5), 3.0)
    totals = np.full((1, 1, 5), 2.0)
    result = uoc_dust_driver(
        start,
        np.array([[[1.2]]]),
        np.array([[[10.0]]]),
        np.array([[[1.5]]]),
        single_cell["ustar"],
        np.array([[1]]),
        single_cell["cf"],
        60.0,
        totals,
        soil,
        single_cell["erod"],
        psd=psd,
    )
    rate = dust_emission_shao2011(psd=psd, **{**single_cell, "gwet": np.array([[1.0]])})
    assert np.array_equal(result["smois"][:, 0, :], [[1.0]])
    assert np.allclose(result["dust"][0, 0, 0], 3.0 + rate[0, 0] * 60.0 / 10.0 / 1.2 * 1e9)
    assert np.allclose(result["emis_dust"], totals + rate * 60.0)
    assert np.isfinite(result["ust_t_min"]).all()
    assert np.isfinite(result["rough_cor"]).all()
    assert np.isfinite(result["smois_cor"]).all()


def test_driver_resets_inactive_totals_and_marks_undefined_diagnostics(single_cell, psd):
    """U:252-266: inactive totals reset; undefined carryover becomes NaN."""
    soil = np.moveaxis(single_cell["soil_frac"], -1, 1)
    dust = np.ones((1, 1, 1, 5))
    totals = np.ones((1, 1, 5))
    kwargs = dict(
        dust=dust,
        rho_phy=np.ones((1, 1, 1)),
        dz8w=np.ones((1, 1, 1)),
        smois=np.full((1, 1, 1), 0.02),
        ust=single_cell["ustar"],
        isltyp=np.array([[2]]),
        vegfra=single_cell["cf"],
        dt=1.0,
        emis_dust=totals,
        soil_top_cat=soil,
        erod=single_cell["erod"],
        psd=psd,
    )
    with pytest.warns(ShaoDustWarning, match="smois_cor"):
        missing = uoc_dust_driver(**kwargs)
    assert np.isnan(missing["smois_cor"][0, 0])
    inactive = uoc_dust_driver(**{**kwargs, "erod": np.zeros((1, 1, 3))})
    assert np.array_equal(inactive["dust"], dust)
    assert np.array_equal(inactive["emis_dust"], np.zeros_like(totals))
    assert inactive["ust_t_min"][0, 0] == -9999.0
    assert inactive["rough_cor"][0, 0] == inactive["smois_cor"][0, 0] == 1.0
    with pytest.warns(ShaoDustWarning, match="all 12"):
        blank = uoc_dust_driver(**{**kwargs, "soil_top_cat": np.zeros_like(soil)})
    assert np.isnan(blank["ust_t_min"][0, 0])
    assert np.isnan(blank["rough_cor"][0, 0])
    assert np.isnan(blank["smois_cor"][0, 0])


def test_default_coefficient_and_all_surface_options(single_cell, psd):
    """S2/S5: source cys is default; host branches share their stated tables."""
    base = dust_emission_shao2011(psd=psd, **single_cell)
    assert np.array_equal(
        base, dust_emission_shao2011(psd=psd, settings=Shao2011Settings(c_y=1.0e-5), **single_cell)
    )
    for sf in (-1, 2, 4):
        assert np.array_equal(
            base, dust_emission_shao2011(psd=psd, sf_surface_physics=sf, **single_cell)
        )
    fallback = dust_emission_shao2011(psd=psd, sf_surface_physics=0, **single_cell)
    for sf in (1, 5, 7, 8):
        assert np.array_equal(
            fallback, dust_emission_shao2011(psd=psd, sf_surface_physics=sf, **single_cell)
        )
    with pytest.warns(ShaoDustWarning, match="unknown sf_surface_physics"):
        unknown = dust_emission_shao2011(psd=psd, sf_surface_physics=99, **single_cell)
    assert np.array_equal(unknown, base)


def test_zero_or_negative_soil_sum_keeps_wrf_arithmetic(single_cell, psd):
    """S4: active zero and negative totals warn without a positive-sum gate."""
    zero = np.zeros((1, 1, 12))
    zero[0, 0, :2] = [0.5, -0.5]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ShaoDustWarning)
        with np.errstate(all="ignore"):
            result = dust_emission_shao2011(psd=psd, **{**single_cell, "soil_frac": zero})
    assert any("sum is not positive" in str(item.message) for item in caught)
    assert np.isnan(result).any() or np.isinf(result).any()
    negative = np.zeros((1, 1, 12))
    negative[0, 0, 0] = -1.0
    with pytest.warns(ShaoDustWarning, match="sum is not positive"):
        dust_emission_shao2011(psd=psd, **{**single_cell, "soil_frac": negative})


def test_zero_soil_sum_below_wind_threshold_still_divides_by_zero(single_cell, psd):
    """Q:188-195,349: a nonzero class with zero total still adds 0/0."""
    zero = np.zeros((1, 1, 12))
    zero[0, 0, :2] = [0.5, -0.5]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ShaoDustWarning)
        with np.errstate(all="ignore"):
            result = dust_emission_shao2011(
                psd=psd, **{**single_cell, "soil_frac": zero, "ustar": np.array([[0.01]])}
            )
    assert np.isnan(result).all()


def test_threshold_minimum_starts_at_999(single_cell, psd):
    """S13: the reported last-class threshold respects Q:173's 999 sentinel."""
    reference = threshold_velocity_shao_lu(psd["dm"], RHOP_DEFAULT, 1.0).min()
    rho = (reference / 1000.5) ** 2
    flux, _, threshold, _, _ = _flux_and_diagnostics(
        np.array([[1000.51]]),
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        single_cell["erod"],
        single_cell["soil_frac"],
        np.array([[rho]]),
        psd,
        2,
        3,
    )
    assert threshold[0, 0] == 999.0
    assert flux.sum() > 0.0


def test_diagnostics_use_last_class_and_selected_isltyp(single_cell, psd):
    """Q:173-185: threshold follows last nonzero class; moisture follows isltyp."""
    soil = np.zeros((1, 1, 12))
    soil[0, 0, :2] = [0.5, 0.5]
    args = (
        single_cell["ustar"],
        single_cell["gwet"],
        single_cell["cf"],
        single_cell["erod"],
        soil,
        single_cell["rhoa"],
        psd,
        2,
        3,
    )
    _, _, threshold, rough, moist = _flux_and_diagnostics(*args, isltyp=np.array([[1]]))
    last = soil.copy()
    last[0, 0, :2] = [0.0, 1.0]
    _, _, last_threshold, _, _ = _flux_and_diagnostics(
        *args[:4], last, *args[5:], isltyp=np.array([[2]])
    )
    assert np.array_equal(threshold, last_threshold)
    assert rough[0, 0] > 1.0
    assert (
        moist[0, 0]
        == moisture_correction_fecan(single_cell["gwet"], NOAH_DRYSMC[0], CC2ISL[0])[0, 0]
    )


def test_inactive_invalid_values_do_not_warn_or_change_moisture(single_cell, psd):
    """Warnings and the input moisture clamp apply only where WRF runs."""
    data = {
        **single_cell,
        "gwet": np.array([[2.0]]),
        "cf": np.array([[200.0]]),
        "erod": np.zeros((1, 1, 3)),
        "soil_frac": -single_cell["soil_frac"],
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ShaoDustWarning)
        _, moisture, *_ = _flux_and_diagnostics(
            data["ustar"],
            data["gwet"],
            data["cf"],
            data["erod"],
            data["soil_frac"],
            data["rhoa"],
            psd,
            99,
            3,
        )
    assert moisture[0, 0] == 2.0
    assert not any(isinstance(item.message, ShaoDustWarning) for item in caught)


def test_active_nonfinite_cover_uses_one_controllable_warning_class(single_cell, psd):
    """An impossible active cover must be filterable as ShaoDustWarning."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        dust_emission_shao2011(psd=psd, **{**single_cell, "cf": np.array([[np.nan]])})
    assert caught and all(isinstance(item.message, ShaoDustWarning) for item in caught)


def test_driver_accepts_twelve_and_sixteen_soil_slots(single_cell, psd):
    """U:209-212: the driver ignores slots 13-16 when present."""
    soil12 = np.moveaxis(single_cell["soil_frac"], -1, 1)
    soil16 = np.concatenate([soil12, np.ones((1, 4, 1))], axis=1)
    inputs = dict(
        dust=np.zeros((1, 1, 1, 5)),
        rho_phy=np.ones((1, 1, 1)),
        dz8w=np.ones((1, 1, 1)),
        smois=np.array([[[0.02]]]),
        ust=single_cell["ustar"],
        isltyp=np.array([[1]]),
        vegfra=single_cell["cf"],
        dt=1.0,
        emis_dust=np.zeros((1, 1, 5)),
        erod=single_cell["erod"],
        psd=psd,
    )
    a = uoc_dust_driver(**inputs, soil_top_cat=soil12)
    b = uoc_dust_driver(**inputs, soil_top_cat=soil16)
    assert all(np.array_equal(a[key], b[key]) for key in a)


def test_driver_rejects_wrong_array_shape_with_value_error(single_cell, psd):
    """Malformed WRF arrays fail at the driver boundary with a clear error."""
    with pytest.raises(ValueError, match="shapes"):
        uoc_dust_driver(
            np.array(0.0),
            np.ones((1, 1, 1)),
            np.ones((1, 1, 1)),
            np.ones((1, 1, 1)),
            single_cell["ustar"],
            np.array([[1]]),
            single_cell["cf"],
            1.0,
            np.zeros((1, 1, 5)),
            np.zeros((1, 12, 1)),
            single_cell["erod"],
            psd=psd,
        )
