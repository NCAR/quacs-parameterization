"""Tests for quacs.dust.ginoux."""

import warnings
from dataclasses import replace

import numpy as np
import pytest

from quacs.dust.ginoux.ginoux import (
    CH_DUST,
    DEN_DUST,
    FRAC_S,
    GRAV,
    IPOINT,
    LAND,
    LAND_MAX,
    NDCLS,
    REFF_DUST,
    GinouxSettings,
    compute_gwet,
    dust_emission,
    ginoux_dust_driver,
    initialize,
    threshold_velocity,
    tunable_names,
)


def _dust_emission_scalar(w10m, u_t, gwet, oro, erod, ch_dust, month):
    """Slow scalar transcription of source_du for the vectorized test."""
    nx, ny = w10m.shape
    n_bins = u_t.shape[-1]
    F = np.zeros((nx, ny, n_bins))
    for n in range(n_bins):
        m = IPOINT[n]
        for i in range(nx):
            for j in range(ny):
                if not (oro[i, j] < LAND_MAX):
                    continue
                if gwet[i, j] < 0.5:
                    u_ts = max(0.0, u_t[i, j, n] * (1.2 + 0.2 * np.log10(max(1.0e-3, gwet[i, j]))))
                else:
                    u_ts = 100.0  # wet surface, no erosion
                srce = FRAC_S[n] * erod[i, j, m]
                dsrc = ch_dust[n, month - 1] * srce * w10m[i, j] ** 2 * (w10m[i, j] - u_ts)
                if dsrc < 0.0:
                    dsrc = 0.0
                F[i, j, n] = dsrc
    return F


@pytest.fixture
def cell():
    """Single land cell, wind above threshold, dry enough to erode."""
    rhoa = np.array([[1.2]])
    return dict(
        w10m=np.array([[12.0]]),
        u_t=threshold_velocity(REFF_DUST, DEN_DUST, rhoa),
        gwet=np.array([[0.1]]),
        oro=np.array([[LAND]]),
        erod=np.full((1, 1, NDCLS), 0.5),
        settings=GinouxSettings(ch_dust=CH_DUST.copy()),
        month=1,
    )


def test_explicit_default_settings_equal_implicit_defaults(cell):
    case = {key: value for key, value in cell.items() if key != "settings"}
    np.testing.assert_array_equal(
        dust_emission(**case), dust_emission(**case, settings=GinouxSettings())
    )


def _scale(value):
    return value * 1.2


def _scale_first(value):
    return (value[0] * 1.2, *value[1:])


def _scale_one_ch_dust(value):
    table = np.asarray(value).copy()
    table[0, 0] *= 1.2
    return tuple(map(tuple, table))


# One entry per setting. A setting without an entry fails the test below.
GINOUX_SETTING_CHANGES = {
    "reff_dust": _scale_first,
    "den_dust": _scale_first,
    "frac_s": _scale_first,
    "ipoint": lambda value: (1, *value[1:]),
    "porosity": lambda value: (*value[:2], value[2] * 1.2, *value[3:]),
    "ch_dust": _scale_one_ch_dust,
    "layer_cutoff": lambda value: 18.0,
    "wet_cutoff": lambda value: 0.3,
    "wet_threshold": _scale,
    "moisture_base": _scale,
    "moisture_log": _scale,
    "moisture_floor": _scale,
    "threshold_scale": _scale,
    "threshold_cohesion": _scale,
    "threshold_drag": _scale,
    "threshold_size": _scale,
    "threshold_offset": _scale,
    "threshold_size_power": _scale,
    "threshold_drag_power": _scale,
    "threshold_diameter_power": _scale,
}


def _ginoux_driver(settings, wind=20.0, dz=15.0, smois=0.02):
    return ginoux_dust_driver(
        dust=np.zeros((1, 1, 1, 5)),
        rho_phy=np.array([[[1.2]]]),
        dz8w=np.array([[[dz]]]),
        dt=1.0,
        u10=np.full((1, 1), wind),
        v10=np.zeros((1, 1)),
        u_phy=np.array([[[20.0]]]),
        v_phy=np.zeros((1, 1, 1)),
        smois=np.array([[[smois]]]),
        isltyp=np.array([[3]]),
        oro=np.ones((1, 1)),
        erod=np.array([[[0.2, 0.7, 0.4]]]),
        emis_dust=np.zeros((1, 1, 5)),
        dustin=np.zeros((1, 1, 5)),
        settings=settings,
    )


def _ginoux_outputs(settings):
    erod = np.array([[[0.2, 0.7, 0.4]]])
    ut = threshold_velocity(settings.reff_dust, settings.den_dust, np.array([[1.2]]), settings)
    pieces = [ut.ravel()]
    pieces.append(compute_gwet(np.array([[0.02]]), np.array([[3]]), settings).ravel())
    for wet, wind in ((0.1, 20.0), (0.0001, 20.0), (0.4, 101.0), (0.6, 101.0)):
        pieces.append(
            dust_emission(
                np.array([[wind]]), ut, np.array([[wet]]), np.ones((1, 1)), erod, settings=settings
            ).ravel()
        )
    pieces.append(_ginoux_driver(settings)["flux"].ravel())
    pieces.append(_ginoux_driver(settings, wind=5.0)["flux"].ravel())
    return np.concatenate(pieces)


@pytest.mark.parametrize("name", tunable_names(GinouxSettings))
def test_each_setting_reaches_an_output_so_tuning_it_is_never_silently_ignored(name):
    base = GinouxSettings()
    changed = GINOUX_SETTING_CHANGES[name](getattr(base, name))
    assert not np.array_equal(
        _ginoux_outputs(replace(base, **{name: changed})), _ginoux_outputs(base)
    ), name


def test_every_setting_case_names_a_real_setting():
    assert set(GINOUX_SETTING_CHANGES) <= set(tunable_names(GinouxSettings))


def test_driver_applies_tuned_porosity_and_cutoff_like_the_standalone_functions():
    """The driver must read porosity and the thin-layer cutoff from settings, as the pieces do."""
    base = GinouxSettings()
    tuned = replace(
        base,
        porosity=(*base.porosity[:2], base.porosity[2] * 0.5, *base.porosity[3:]),
        layer_cutoff=18.0,
    )
    erod = np.array([[[0.2, 0.7, 0.4]]])
    for dz, wind in ((15.0, 20.0), (10.0, 20.0)):
        thin = dz < tuned.layer_cutoff
        gwet = compute_gwet(np.array([[0.02]]), np.array([[3]]), tuned)
        ut = threshold_velocity(tuned.reff_dust, tuned.den_dust, np.array([[1.2]]), tuned)
        expected = dust_emission(
            np.array([[20.0 if thin else wind]]), ut, gwet, np.ones((1, 1)), erod, settings=tuned
        )
        np.testing.assert_array_equal(_ginoux_driver(tuned, wind=wind, dz=dz)["flux"], expected)
    assert tuned.layer_cutoff != base.layer_cutoff


def test_default_numbers_are_pinned_so_a_default_drift_fails(cell):
    """Defaults are the WRF values; a 1e-6 drift in any of them must show here."""
    rhoa = np.array([[1.2]])
    np.testing.assert_allclose(
        threshold_velocity(REFF_DUST, DEN_DUST, rhoa).ravel(),
        [
            2.503409201494799,
            1.5354480450676777,
            1.0242925179416447,
            0.6393881267105276,
            0.4185922464702874,
        ],
        rtol=1e-9,
        atol=0,
    )
    flux = dust_emission(**{k: v for k, v in cell.items() if k != "settings"})
    assert flux.min() > 0
    np.testing.assert_allclose(
        flux.ravel(),
        [
            5.4700362999389964e-08,
            1.5068954815102543e-07,
            1.5805018774164032e-07,
            1.6359281097536842e-07,
            1.6677227165082786e-07,
        ],
        rtol=1e-9,
        atol=0,
    )
    driver = _ginoux_driver(None)["flux"].ravel()
    assert driver.min() > 0
    np.testing.assert_allclose(
        driver,
        [
            2.2611264098860778e-07,
            1.0398010136107253e-06,
            1.0664995107005286e-06,
            1.0866037025243267e-06,
            1.0981362358790548e-06,
        ],
        rtol=1e-9,
        atol=0,
    )


# --- shape / contract ---------------------------------------------------
def test_output_shape(cell):
    assert dust_emission(**cell).shape == (1, 1, 5)


def test_all_fluxes_nonnegative(cell):
    assert (dust_emission(**cell) >= 0.0).all()


def test_signature_has_no_frlake():
    import inspect

    params = inspect.signature(dust_emission).parameters
    assert "frlake" not in params, "dust_opt=1 has no lake-fraction term"


# --- the core fix: per-bin threshold -----------------------------------
def test_per_bin_threshold_is_used(cell):
    """
    Flux ratio between bins must follow frac_s_i*(w-u_t_i) / frac_s_j*(w-u_t_j),
    i.e. the per-bin threshold must enter the flux (not a single bin mean).
    """
    c = cell
    F = dust_emission(**c)[0, 0, :]
    w = c["w10m"][0, 0]
    ut = c["u_t"][0, 0, :] * (1.2 + 0.2 * np.log10(c["gwet"][0, 0]))
    ut = np.maximum(ut, 0.0)
    expect = FRAC_S * (w - ut)  # erod equal across bins here
    expect = expect / expect[0]
    got = F / F[0]
    assert np.allclose(got, expect, rtol=1e-12)


def test_bin_mean_would_differ(cell):
    """Guard: a single bin-mean threshold gives a different size split."""
    c = cell
    F_binwise = dust_emission(**c)[0, 0, :]
    c_mean = dict(c)
    c_mean["u_t"] = np.full_like(c["u_t"], c["u_t"].mean())
    F_mean = dust_emission(**c_mean)[0, 0, :]
    split_bw = F_binwise / F_binwise.sum()
    split_mn = F_mean / F_mean.sum()
    assert not np.allclose(split_bw, split_mn, atol=1e-4)


# --- vectorized vs scalar transcription of source_du -------------------
def test_vectorized_matches_scalar():
    rng = np.random.default_rng(0)
    nx, ny = 9, 7
    rhoa = rng.uniform(0.9, 1.3, (nx, ny))
    kw = dict(
        w10m=rng.uniform(0.0, 25.0, (nx, ny)),
        u_t=threshold_velocity(REFF_DUST, DEN_DUST, rhoa),
        gwet=rng.uniform(0.0, 0.7, (nx, ny)),
        oro=rng.choice([LAND, 2.0], (nx, ny)),
        erod=rng.uniform(0.0, 1.0, (nx, ny, NDCLS)),
        settings=GinouxSettings(),
        month=1,
    )
    F_vec = dust_emission(**kw)
    F_scl = _dust_emission_scalar(
        **{key: value for key, value in {**kw, "ch_dust": CH_DUST}.items() if key != "settings"}
    )
    assert np.allclose(F_vec, F_scl, rtol=0.0, atol=1e-300) or np.max(
        np.abs(F_vec - F_scl)
    ) < 1e-13 * max(1.0, np.max(np.abs(F_scl)))


# --- threshold formula vs WRF CGS expression ---------------------------
def test_threshold_matches_wrf_cgs():
    def u_ts0_wrf(reff_m, den_kgm3, rhoa_kgm3, g=GRAV):
        g0 = g * 1.0e2
        den = den_kgm3 * 1.0e-3
        diam = 2.0 * reff_m * 1.0e2
        rhoa = rhoa_kgm3 * 1.0e-3
        return (
            0.13
            * 1.0e-2
            * np.sqrt(den * g0 * diam / rhoa)
            * np.sqrt(1.0 + 0.006 / den / g0 / diam**2.5)
            / np.sqrt(1.928 * (1331.0 * diam**1.56 + 0.38) ** 0.092 - 1.0)
        )

    rhoa = 1.2
    wrf = np.array([u_ts0_wrf(r, d, rhoa) for r, d in zip(REFF_DUST, DEN_DUST)])
    py = threshold_velocity(REFF_DUST, DEN_DUST, np.array([[rhoa]]))[0, 0, :]
    assert np.max(np.abs(wrf - py)) < 1e-13


def test_threshold_decreases_with_size():
    u_t = threshold_velocity(REFF_DUST, DEN_DUST, np.full((1, 1), 1.2))
    assert np.all(np.diff(u_t[0, 0, :]) < 0)


# --- cutoffs -----------------------------------------------------------
def test_zero_over_ocean(cell):
    c = dict(cell)
    c["oro"] = np.array([[2.0]])
    assert (dust_emission(**c) == 0.0).all()


def test_zero_wet_soil(cell):
    c = dict(cell)
    c["gwet"] = np.array([[0.5]])
    assert (dust_emission(**c) == 0.0).all()


def test_wet_wind_above_100_emits_because_wrf_uses_a_threshold(cell):
    c = {**cell, "gwet": np.array([[0.5]]), "w10m": np.array([[101.0]])}
    assert (dust_emission(**c) > 0).all()


def test_negative_coefficient_below_threshold_emits_after_product_clipping(cell):
    c = {**cell, "w10m": np.array([[0.2]]), "settings": replace(cell["settings"], ch_dust=-CH_DUST)}
    with pytest.warns(UserWarning, match="ch_dust"):
        assert (dust_emission(**c) > 0).all()


def test_ch_dust_selects_each_bin_at_start_month(cell):
    table = CH_DUST.copy()
    table[:, 1] = table[:, 0] * np.arange(1, 6)
    base = dust_emission(**cell)
    selected = dust_emission(
        **{**cell, "settings": replace(cell["settings"], ch_dust=table), "month": 2}
    )
    np.testing.assert_allclose(selected[0, 0], base[0, 0] * np.arange(1, 6))


def test_sea_skips_invalid_inputs_without_warning(cell):
    c = {
        **cell,
        "oro": np.array([[2.0]]),
        "gwet": np.array([[np.nan]]),
        "erod": np.full((1, 1, 3), np.nan),
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert (dust_emission(**c) == 0).all()


def test_land_erod_warning_keeps_wrf_arithmetic(cell):
    c = {**cell, "erod": np.full((1, 1, 3), -0.2)}
    with pytest.warns(UserWarning, match="erod"):
        assert (dust_emission(**c) >= 0).all()


def test_land_gwet_warning_keeps_wrf_threshold(cell):
    """Saturation above one warns while source_du still uses its wet threshold."""
    c = {**cell, "gwet": np.array([[1.2]]), "w10m": np.array([[101.0]])}
    with pytest.warns(UserWarning, match="gwet"):
        assert (dust_emission(**c) > 0).all()


def test_wrong_ch_dust_shape_is_an_error(cell):
    """A missing monthly column cannot be interpreted as WRF ch_dust."""
    with pytest.raises(ValueError, match="ch_dust"):
        dust_emission(**{**cell, "settings": replace(cell["settings"], ch_dust=np.ones(5))})


def test_initialize_warns_on_invalid_land_erod(tmp_path):
    """Static land erodibility warnings retain the loaded WRF values."""
    e_f = tmp_path / "erod.npy"
    o_f = tmp_path / "oro.npy"
    np.save(e_f, np.full((1, 1, 3), -0.2))
    np.save(o_f, np.ones((1, 1)))
    with pytest.warns(UserWarning, match="erod"):
        erod, _ = initialize(e_f, o_f)
    assert erod[0, 0, 0] == -0.2


def test_driver_rejects_invalid_land_soil_type(cell):
    """The WRF porosity table has only soil types 1 through 19."""
    args = dict(
        dust=np.zeros((1, 2, 1, 5)),
        rho_phy=np.ones((1, 2, 1)),
        dz8w=np.ones((1, 2, 1)),
        dt=1.0,
        u10=np.ones((1, 1)),
        v10=np.zeros((1, 1)),
        u_phy=np.ones((1, 2, 1)),
        v_phy=np.zeros((1, 2, 1)),
        smois=np.ones((1, 1, 1)),
        isltyp=np.array([[20]]),
        oro=np.ones((1, 1)),
        erod=np.ones((1, 1, 3)),
        emis_dust=np.zeros((1, 1, 5)),
        dustin=np.zeros((1, 1, 5)),
    )
    with pytest.raises(ValueError, match="isltyp"):
        ginoux_dust_driver(**args)
    args["oro"][:] = 2.0
    assert (ginoux_dust_driver(**args)["flux"] == 0).all()


def test_driver_thin_layer_uses_model_wind_and_updates_accumulators(cell):
    args = dict(
        dust=np.ones((1, 2, 1, 5)),
        rho_phy=np.ones((1, 2, 1)),
        dz8w=np.full((1, 2, 1), 11.0),
        dt=2.0,
        u10=np.array([[0.0]]),
        v10=np.array([[0.0]]),
        u_phy=np.full((1, 2, 1), 20.0),
        v_phy=np.zeros((1, 2, 1)),
        smois=np.full((1, 1, 1), 0.01),
        isltyp=np.array([[3]]),
        oro=np.array([[1.0]]),
        erod=np.full((1, 1, 3), 0.5),
        emis_dust=np.ones((1, 1, 5)),
        dustin=np.ones((1, 1, 5)),
    )
    out = ginoux_dust_driver(**args)
    assert (out["flux"] > 0).all()
    np.testing.assert_allclose(out["emis_dust"], 1 + out["flux"] * 2)
    np.testing.assert_allclose(out["dust"][:, 0], 1 + out["flux"] * 2 / 11 * 1e9)
    args["dz8w"][:] = 12.0
    assert (ginoux_dust_driver(**args)["flux"] == 0).all()


def test_driver_chem_opt_uses_tiny_start_and_writes_dustin(cell):
    args = dict(
        dust=np.full((1, 2, 1, 5), 3.0),
        rho_phy=np.ones((1, 2, 1)),
        dz8w=np.full((1, 2, 1), 50.0),
        dt=1.0,
        u10=np.array([[20.0]]),
        v10=np.array([[0.0]]),
        u_phy=np.zeros((1, 2, 1)),
        v_phy=np.zeros((1, 2, 1)),
        smois=np.full((1, 1, 1), 0.01),
        isltyp=np.array([[3]]),
        oro=np.array([[1.0]]),
        erod=np.full((1, 1, 3), 0.5),
        emis_dust=np.zeros((1, 1, 5)),
        dustin=np.zeros((1, 1, 5)),
    )
    for chem_opt in (2, 11):
        out = ginoux_dust_driver(**args, chem_opt=chem_opt)
        np.testing.assert_array_equal(out["dust"], args["dust"])
        np.testing.assert_allclose(out["dustin"], 1e-16 + out["flux"] / 50 * 1e9)


def test_negative_air_density_gives_zero_threshold_like_gfortran(cell):
    """
    Negative air density makes sqrt(.../rhoa) NaN. gfortran MAX(0, NaN) is 0
    (F L193), so WRF treats the threshold as 0 and the land cell still
    emits a finite flux. np.maximum would return NaN and poison the tracers.
    A comparison with the compiled pinned Fortran covers the same cells.
    """
    args = dict(
        dust=np.zeros((1, 2, 1, 5)),
        rho_phy=np.full((1, 2, 1), -1.0),
        dz8w=np.full((1, 2, 1), 50.0),
        dt=2.0,
        u10=np.array([[12.0]]),
        v10=np.zeros((1, 1)),
        u_phy=np.zeros((1, 2, 1)),
        v_phy=np.zeros((1, 2, 1)),
        smois=np.full((1, 1, 1), 0.03),
        isltyp=np.array([[3]]),
        oro=np.ones((1, 1)),
        erod=np.full((1, 1, 3), 0.3),
        emis_dust=np.zeros((1, 1, 5)),
        dustin=np.zeros((1, 1, 5)),
    )
    with np.errstate(invalid="ignore"):
        out = ginoux_dust_driver(**args)
    expect = CH_DUST[:, 0] * FRAC_S * 0.3 * 12.0**3
    assert np.all(np.isfinite(out["flux"])) and np.all(np.isfinite(out["dust"]))
    np.testing.assert_allclose(out["flux"][0, 0], expect, rtol=1e-12)


def test_zero_below_threshold(cell):
    c = dict(cell)
    c["w10m"] = np.array([[0.2]])
    assert (dust_emission(**c) == 0.0).all()


def test_moisture_raises_threshold(cell):
    dry = dict(cell)
    dry["gwet"] = np.array([[0.01]])
    moist = dict(cell)
    moist["gwet"] = np.array([[0.40]])
    assert dust_emission(**moist).sum() < dust_emission(**dry).sum()


# --- scaling -----------------------------------------------------------
def test_scales_linearly_with_ch_dust(cell):
    F1 = dust_emission(**cell)
    F2 = dust_emission(
        **{**cell, "settings": replace(cell["settings"], ch_dust=cell["settings"].ch_dust * 2.0)}
    )
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-12)


def test_scales_linearly_with_erod(cell):
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, "erod": cell["erod"] * 2.0})
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-12)


# --- source-verified parameters ----------------------------------------
def test_params_match_source():
    assert np.allclose(REFF_DUST, np.array([0.73, 1.4, 2.4, 4.5, 8.0]) * 1e-6)
    assert DEN_DUST[0] == 2500.0 and np.all(DEN_DUST[1:] == 2650.0)
    assert np.allclose(FRAC_S, [0.1, 0.25, 0.25, 0.25, 0.25])
    assert np.isclose(FRAC_S.sum(), 1.1)
    assert np.array_equal(IPOINT, np.array([3, 2, 2, 2, 2]) - 1)
    assert CH_DUST.shape == (5, 12) and np.all(CH_DUST == 0.8e-9)


def test_sand_class_unused():
    """ipoint never selects class 0 (Sand), per WRF ipoint=(3,2,2,2,2)."""
    assert 0 not in IPOINT


# --- vectorized grid ---------------------------------------------------
def test_vectorized_grid():
    nx, ny, nb = 12, 10, 5
    rng = np.random.default_rng(1)
    rhoa = rng.uniform(0.9, 1.3, (nx, ny))
    F = dust_emission(
        w10m=rng.uniform(0, 20, (nx, ny)),
        u_t=threshold_velocity(REFF_DUST, DEN_DUST, rhoa),
        gwet=rng.uniform(0, 0.45, (nx, ny)),
        oro=rng.choice([LAND, 2.0], (nx, ny)),
        erod=rng.uniform(0, 1, (nx, ny, NDCLS)),
    )
    assert F.shape == (nx, ny, nb)
    assert (F >= 0.0).all()


# --- initialize() ------------------------------------------------------
def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("no_erod.npy", "no_oro.npy")


def test_initialize_valid(tmp_path):
    e_f = tmp_path / "erod.npy"
    o_f = tmp_path / "oro.npy"
    np.save(e_f, np.full((6, 6, NDCLS), 0.4))
    np.save(o_f, np.ones((6, 6)))
    erod, oro = initialize(str(e_f), str(o_f))
    assert erod.shape == (6, 6, NDCLS) and oro.shape == (6, 6)
