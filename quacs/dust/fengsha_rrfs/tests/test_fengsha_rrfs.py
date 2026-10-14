"""Tests for quacs.dust.fengsha_rrfs."""

import warnings
from dataclasses import replace

import numpy as np
import pytest

import quacs.dust.fengsha_rrfs.fengsha_rrfs as rrfs
from quacs.dust.fengsha_rrfs.fengsha_rrfs import (
    ALPHA_DEFAULT,
    CON_RD,
    DRY_FACTOR_DEFAULT,
    DUST12M_FIELDS,
    GAMMA_DEFAULT,
    GRAV,
    KOK_FRAC,
    KVHMAX_DEFAULT,
    MOIST_CORR_DEFAULT,
    N_BINS,
    RHO_DRY,
    RHOSOIL,
    UTHRES_TABLE,
    FengshaRRFSWarning,
    RrfsSettings,
    RrfsStepResult,
    dust_emission,
    fengsha_rrfs_driver,
    fengsha_rrfs_monthly_driver,
    flux_ratio_mb95,
    initialize,
    moisture_correction_fecan,
    moisture_correction_shao,
    rrfs_dust_step,
    select_month,
    soil_moisture_vol2grav,
    tunable_names,
)


@pytest.fixture
def host_case():
    """Host arrays with distinct months and vertical layers for adapter tests."""
    dust12m = np.zeros((1, 1, 12, 5))
    dust12m[0, 0, 2] = [0.1, 0.6, 0.5, 0.8, 0.2]
    dust12m[0, 0, 3] = [0.15, 0.6, 0.55, 0.8, 0.18]
    return dict(
        dust12m=dust12m,
        month=4,
        zorl=np.array([[1.0]]),
        snow=np.array([[0.0]]),
        land=np.array([[1]], dtype=int),
        smc=np.array([[[0.3, 0.02, 0.4]]]),
        tslb=np.array([[[290.0, 250.0, 260.0]]]),
        pr3d=np.array([[[100000.0, 99000.0, 98000.0]]]),
        prl3d=np.array([[[103338.0, 85000.0]]]),
        tk3d=np.array([[[300.0, 270.0]]]),
        ustar=np.array([[1.2]]),
        garea=np.array([[2.0]]),
        soiltyp=np.array([[2]]),
        dt=10.0,
        dust=np.array([[4.0]]),
        coarse_pm=np.array([[7.0]]),
        settings=RrfsSettings(
            alpha=0.02, gamma=1.1, kvhmax=3.0e-4, moist_opt=1, dry_factor=1.4, moist_corr=0.9
        ),
    )


@pytest.fixture
def captured_step(monkeypatch):
    """Observe the adapter arguments while still running the real driver."""
    captured = {}
    driver = fengsha_rrfs_driver

    def capture(**kwargs):
        captured.update(kwargs)
        return driver(**kwargs)

    monkeypatch.setattr(rrfs, "fengsha_rrfs_driver", capture)
    return captured


def test_step_selects_calendar_month_and_five_monthly_fields(host_case, captured_step):
    rrfs_dust_step(**host_case)
    for name, value in zip(DUST12M_FIELDS, host_case["dust12m"][0, 0, 3]):
        np.testing.assert_array_equal(captured_step[name], [[value]])


def test_step_converts_zorl_centimeters_to_znt_meters(host_case, captured_step):
    case = {**host_case, "zorl": np.array([[25.0]])}
    result = rrfs_dust_step(**case)
    np.testing.assert_array_equal(captured_step["znt"], [[0.25]])
    assert not result.rate.any()


def test_step_converts_snow_millimeters_to_snowh_meters(host_case, captured_step):
    case = {**host_case, "snow": np.array([[1.0]])}
    result = rrfs_dust_step(**case)
    np.testing.assert_array_equal(captured_step["snowh"], [[0.001]])
    assert not result.rate.any()


def test_step_converts_land_to_float_xland(host_case, captured_step):
    case = {**host_case, "land": np.array([[2]], dtype=int)}
    result = rrfs_dust_step(**case)
    assert captured_step["xland"].dtype.kind == "f"
    np.testing.assert_array_equal(captured_step["xland"], [[2.0]])
    assert not result.rate.any()


def test_step_uses_soil_moisture_layer_two(host_case, captured_step):
    result = rrfs_dust_step(**host_case)
    np.testing.assert_array_equal(captured_step["slc"], [[0.02]])
    assert result.rate.sum() > 0


def test_step_uses_soil_temperature_layer_one(host_case, captured_step):
    result = rrfs_dust_step(**host_case)
    np.testing.assert_array_equal(captured_step["stemp"], [[290.0]])
    assert result.rate.sum() > 0


def test_step_uses_first_two_pressure_edges_for_airmass(host_case, captured_step):
    rrfs_dust_step(**host_case)
    np.testing.assert_array_equal(captured_step["p8w1"], [[100000.0]])
    np.testing.assert_array_equal(captured_step["p8w2"], [[99000.0]])


def test_step_uses_first_middle_pressure_for_density(host_case, captured_step):
    rrfs_dust_step(**host_case)
    np.testing.assert_array_equal(captured_step["p"], [[103338.0]])


def test_step_uses_first_middle_temperature_for_density(host_case, captured_step):
    rrfs_dust_step(**host_case)
    np.testing.assert_array_equal(captured_step["temperature"], [[300.0]])


def test_step_passes_host_state_and_tuning_without_change(host_case, captured_step):
    result = rrfs_dust_step(**host_case)
    for host_name, driver_name in [
        ("ustar", "ustar"),
        ("garea", "area"),
        ("soiltyp", "isltyp"),
    ]:
        assert captured_step[driver_name] is host_case[host_name]
    for name in ("dt", "settings"):
        assert captured_step[name] == host_case[name]
    assert isinstance(result, RrfsStepResult)


# ===========================================================================
# Pytest test suite   -   run: pytest dust_fengsha_rrfs.py -v
# ===========================================================================


@pytest.fixture
def cell():
    """Single erodible land cell, wind above threshold."""
    return dict(
        ustar=np.array([[0.9]]),
        slc=np.array([[0.05]]),
        rho_a=np.array([[1.2]]),
        clay=np.array([[0.15]]),
        sand=np.array([[0.55]]),
        ssm=np.array([[0.8]]),
        rdrag=np.array([[0.6]]),
        uthrs=np.array([[0.30]]),
        isltyp=np.array([[2]]),
        stemp=np.array([[290.0]]),
        znt=np.array([[0.01]]),
        snowh=np.array([[0.0]]),
        xland=np.array([[1.0]]),
        settings=RrfsSettings(moist_corr=1.0, dry_factor=1.0),
    )


def test_explicit_default_settings_equal_implicit_defaults():
    case = {key: value for key, value in cell.__wrapped__().items() if key != "settings"}
    np.testing.assert_array_equal(
        dust_emission(**case), dust_emission(**case, settings=RrfsSettings())
    )


def _scale(value):
    return value * 1.2


def _scale_first(value):
    return (value[0] * 1.2, *value[1:])


# One entry per setting. A setting without an entry fails the test below.
RRFS_SETTING_CHANGES = {
    "reff_dust": _scale_first,
    "lo_dust": _scale_first,
    "up_dust": _scale_first,
    "mmd_dust": _scale,
    "gsd_dust": _scale,
    "lambda_dust": _scale,
    "cv": _scale,
    "alpha": _scale,
    "gamma": _scale,
    "kvhmax": _scale,
    "moist_opt": lambda value: 2,
    "dry_factor": _scale,
    "moist_corr": _scale,
    "rhow": _scale,
    "rho_dry": _scale,
    "sat_base": _scale,
    "sat_sand": _scale,
    "fecan_clay_a": _scale,
    "fecan_clay_b": _scale,
    "fecan_scale": _scale,
    "fecan_power": _scale,
    "shao_cutoff": lambda value: 0.02,
    "shao_low": _scale,
    "shao_high": _scale,
    "shao_offset": _scale,
    "ratio_cutoff": lambda value: 0.1,
    "ratio_clay": _scale,
    "ratio_offset": _scale,
    "frozen_temp_k": lambda value: 300.0,
    "znt_max": lambda value: 0.005,
    "snow_max": lambda value: -0.001,
    "nonerodible_soiltypes": lambda value: (*value, 2),
    "remap_fine": _scale,
    "remap_coarse": _scale,
    "epsilc": _scale,
    "dust_cap": _scale,
    "coarse_pm_cap": _scale,
}


def _rrfs_cell():
    return {key: value for key, value in cell.__wrapped__().items() if key != "settings"}


def _rrfs_driver_inputs(case):
    driver = {key: value for key, value in case.items() if key != "rho_a"}
    driver.update(
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.ones((1, 1)),
        dt=10.0,
        dust=np.ones((1, 1, 5)),
    )
    return driver


def _host_step_inputs(dust, coarse_pm, ustar=1.2):
    return dict(
        dust12m=np.tile([0.1, 0.6, 0.5, 0.8, 0.2], (1, 1, 12, 1)),
        month=1,
        zorl=np.array([[1.0]]),
        snow=np.array([[0.0]]),
        land=np.array([[1]]),
        smc=np.array([[[0.3, 0.02, 0.4]]]),
        tslb=np.array([[[290.0, 250.0, 260.0]]]),
        pr3d=np.array([[[100000.0, 99000.0, 98000.0]]]),
        prl3d=np.array([[[103338.0, 85000.0]]]),
        tk3d=np.array([[[300.0, 270.0]]]),
        ustar=np.array([[ustar]]),
        garea=np.array([[2.0]]),
        soiltyp=np.array([[2]]),
        dt=10.0,
        dust=np.array([[dust]]),
        coarse_pm=np.array([[coarse_pm]]),
    )


def _rrfs_outputs(settings):
    case = _rrfs_cell()
    case["ustar"] = np.array([[2.0]])
    cases = [
        case,
        {**case, "clay": np.array([[0.3]])},
        {**case, "slc": np.array([[0.025]])},
        {**case, "slc": np.array([[0.04]]), "ustar": np.array([[15.0]])},
    ]
    result = [dust_emission(**sample, settings=settings).ravel() for sample in cases]
    result.append(dust_emission(**case, settings=replace(settings, moist_opt=2)).ravel())
    result.append(dust_emission(**cases[2], settings=replace(settings, moist_opt=2)).ravel())
    result.append(dust_emission(**cases[3], settings=replace(settings, moist_opt=2)).ravel())
    result.append(fengsha_rrfs_driver(**_rrfs_driver_inputs(case), settings=settings).dust.ravel())
    for dust, coarse in ((0.0, 0.0), (300.0, 6000.0)):  # floor, caps (W:589-590,1066-1067)
        inputs = {**_host_step_inputs(dust, coarse), "land": np.array([[2]])}
        step = rrfs_dust_step(**inputs, settings=settings)
        result.append(np.array([step.dust.item(), step.coarse_pm.item()]))
    return np.concatenate(result)


@pytest.mark.parametrize("name", tunable_names(RrfsSettings))
def test_each_setting_reaches_an_output_so_tuning_it_is_never_silently_ignored(name):
    base = RrfsSettings(moist_corr=1.0, dry_factor=1.0)
    changed = RRFS_SETTING_CHANGES[name](getattr(base, name))
    assert not np.array_equal(
        _rrfs_outputs(replace(base, **{name: changed})), _rrfs_outputs(base)
    ), name


def test_every_setting_case_names_a_real_setting():
    assert set(RRFS_SETTING_CHANGES) <= set(tunable_names(RrfsSettings))


@pytest.mark.parametrize(
    "name,excluding_value",
    [
        ("frozen_temp_k", 300.0),
        ("znt_max", 0.005),
        ("snow_max", -0.001),
        ("nonerodible_soiltypes", (2,)),
    ],
)
def test_driver_tracer_gate_reads_the_same_settings_as_the_flux_mask(name, excluding_value):
    """A cell the settings exclude must keep its tracers untouched (no remap), not only get zero flux."""
    case = _rrfs_cell()
    case["ustar"] = np.array([[2.0]])
    inputs = _rrfs_driver_inputs(case)
    default = fengsha_rrfs_driver(**inputs)
    assert default.rate.sum() > 0
    assert not np.array_equal(default.dust, inputs["dust"])
    tuned = fengsha_rrfs_driver(
        **inputs, settings=replace(RrfsSettings(), **{name: excluding_value})
    )
    assert tuned.rate.sum() == 0
    np.testing.assert_array_equal(tuned.dust, inputs["dust"])


def test_default_numbers_are_pinned_so_a_default_drift_fails():
    """Defaults are the RRFS values; a 1e-6 drift (for example sat_base) must show here."""
    base = _rrfs_cell()
    base["ustar"] = np.array([[2.0]])
    fecan = [
        (
            base,
            [
                1.719334911494301e-05,
                1.620854453646707e-05,
                3.32582744070706e-05,
                7.710350223518683e-05,
                1.6316506376709428e-05,
            ],
        ),
        (
            {**base, "clay": np.array([[0.3]])},
            [
                3.5322997807409436e-05,
                3.329975906935583e-05,
                6.832769730349711e-05,
                0.0001584058360119049,
                3.3521562036344405e-05,
            ],
        ),
        (
            {**base, "slc": np.array([[0.025]])},
            [
                1.9633075320658656e-05,
                1.850852754720924e-05,
                3.7977604136569966e-05,
                8.80444442062976e-05,
                1.863180910957762e-05,
            ],
        ),
        (
            {**base, "slc": np.array([[0.04]]), "ustar": np.array([[15.0]])},
            [
                0.007856271312752048,
                0.007406278009712377,
                0.015196924426365341,
                0.03523141586164549,
                0.007455609731106367,
            ],
        ),
    ]
    for sample, expected in fecan:
        flux = dust_emission(**sample).ravel()
        assert flux.min() > 0
        np.testing.assert_allclose(flux, expected, rtol=1e-9, atol=0)
    shao_moisture = [
        (
            0.005,
            [
                2.0594479919219538e-05,
                1.9414864593538098e-05,
                3.983731498995848e-05,
                9.235585910972862e-05,
                1.9544183083847868e-05,
            ],
        ),
        (
            0.01,
            [
                2.047741279927288e-05,
                1.9304502870832066e-05,
                3.9610864030741494e-05,
                9.183087209968594e-05,
                1.9433086264005333e-05,
            ],
        ),
    ]
    for slc, expected in shao_moisture:
        flux = dust_emission(
            **{**base, "slc": np.array([[slc]])}, settings=RrfsSettings(moist_opt=2)
        ).ravel()
        assert flux.min() > 0
        np.testing.assert_allclose(flux, expected, rtol=1e-9, atol=0)
    driver = fengsha_rrfs_driver(**_rrfs_driver_inputs({**base, "slc": np.array([[0.05]])}))
    np.testing.assert_allclose(
        driver.dust.ravel(),
        [
            2141.9789275571966,
            1590.5152327854478,
            3262.5225671409885,
            7562.270601946948,
            11960.421045296745,
        ],
        rtol=1e-9,
        atol=0,
    )


# --- shape and basic contract ---


def test_output_shape(cell):
    assert dust_emission(**cell).shape == (1, 1, N_BINS)


def test_all_fluxes_nonnegative(cell):
    assert (dust_emission(**cell) >= 0.0).all()


def test_bins_sum_to_total(cell):
    """The five rates retain the Kok partition."""
    F_p = dust_emission(**cell)
    np.testing.assert_allclose(F_p[0, 0] / F_p.sum(), KOK_FRAC)


# --- emission cutoffs ---


def test_zero_flux_below_threshold(cell):
    """rustar <= u_thresh must give zero flux (q = 0)."""
    c = {**cell, "ustar": np.array([[0.01]])}
    assert (dust_emission(**c) == 0.0).all()


def test_zero_flux_nonerodible_soiltype(cell):
    """isltyp in {0,15,16,18} must zero the flux."""
    F = dust_emission(**{**cell, "isltyp": np.array([[15]])})
    assert (F == 0.0).all()


def test_soil_type_18_blocks_flux_even_with_strong_wind(cell):
    windy = {**cell, "ustar": np.array([[2.0]])}
    assert dust_emission(**{**windy, "isltyp": np.array([[17]])}).sum() > 0
    assert not dust_emission(**{**windy, "isltyp": np.array([[18]])}).any()


def test_zero_flux_frozen_soil(cell):
    """stemp below 268 K must zero the flux."""
    F = dust_emission(**{**cell, "stemp": np.array([[260.0]])})
    assert (F == 0.0).all()


def test_physical_and_namelist_constants_are_pinned():
    """Defaults are the real RRFS forecast namelist values (input.nml:190-194)."""
    assert (GRAV, CON_RD) == (9.80665, 287.05)
    assert (ALPHA_DEFAULT, GAMMA_DEFAULT, MOIST_CORR_DEFAULT, DRY_FACTOR_DEFAULT) == (
        10.0,
        1.3,
        2.0,
        0.5,
    )


@pytest.mark.parametrize(
    "gate,value",
    [
        ("znt", 0.201),
        ("snowh", 0.001),
        ("xland", 1.5),
        ("stemp", 267.99),
        ("isltyp", 0),
    ],
)
def test_driver_gates_suppress_emission(gate, value, cell):
    assert not dust_emission(**{**cell, gate: np.array([[value]])}).any()


def test_gate_boundaries_still_emit(cell):
    for gate, value in [("znt", 0.2), ("snowh", 0.0), ("xland", 1.49), ("stemp", 268.0)]:
        assert dust_emission(**{**cell, gate: np.array([[value]])}).sum() > 0


def test_soil_fields_and_threshold_are_read_each_call(cell):
    high = dust_emission(**{**cell, "uthrs": np.array([[2.0]])})
    low = dust_emission(**{**cell, "ssm": np.array([[0.4]])})
    assert high.sum() == 0
    assert low.sum() < dust_emission(**cell).sum()


def test_out_of_range_soil_warns_and_keeps_fortran_calculation(cell):
    bad = {**cell, "clay": np.array([[1.2]])}
    with pytest.warns(FengshaRRFSWarning):
        actual = dust_emission(**bad)
    assert np.isfinite(actual).all()
    assert actual.sum() > 0


def test_inactive_invalid_soil_does_not_warn(cell):
    bad = {
        **cell,
        "clay": np.array([[0.15, 1.2]]),
        "sand": np.array([[0.55, 0.3]]),
        "ssm": np.array([[0.8, 1.2]]),
        "settings": replace(cell["settings"], alpha=np.array([[0.01, -0.01]])),
        "xland": np.array([[1.0, 2.0]]),
    }
    for key, value in cell.items():
        if key not in bad:
            bad[key] = value if key == "settings" else np.repeat(value, 2, axis=1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", FengshaRRFSWarning)
        flux = dust_emission(**bad)
    assert caught == []
    assert flux[0, 0].sum() > 0
    assert not flux[0, 1].any()


def test_shao_switches_branch_at_exactly_point_zero_three():
    actual = moisture_correction_shao(np.array([0.029999, 0.03, 0.030001]))
    expected = np.array(
        [np.exp(22.7 * 0.029999), np.exp(95.3 * 0.03 - 2.029), np.exp(95.3 * 0.030001 - 2.029)]
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


def test_fecan_uses_point_six_eight_exponent():
    grav = soil_moisture_vol2grav(0.08, 0.5)
    dry = 0.1 * (14.0 * 0.1 + 17.0)
    expected = np.sqrt(1.0 + 1.21 * max(0.0, grav - dry) ** 0.68)
    np.testing.assert_allclose(
        moisture_correction_fecan(0.08, 0.5, 0.1, RrfsSettings(dry_factor=1.0)),
        expected,
        rtol=1e-13,
    )


def test_moist_correction_changes_layer_two_threshold(cell):
    windy = {**cell, "ustar": np.array([[2.0]])}
    adjusted = dust_emission(**{**windy, "settings": replace(cell["settings"], moist_corr=1.4)})
    expected = dust_emission(**{**windy, "slc": windy["slc"] * 1.4})
    np.testing.assert_allclose(adjusted, expected, rtol=1e-13)
    assert adjusted.sum() < dust_emission(**windy).sum()


def test_drylimit_factor_changes_fecan_threshold(cell):
    windy = {**cell, "ustar": np.array([[2.0]])}
    adjusted = dust_emission(**{**windy, "settings": replace(cell["settings"], dry_factor=2.0)})
    assert adjusted.sum() > dust_emission(**windy).sum()


def test_negative_alpha_clamps_each_bin_after_warning(cell):
    with pytest.warns(FengshaRRFSWarning):
        assert not dust_emission(
            **{**cell, "settings": replace(cell["settings"], alpha=-1.0)}
        ).any()


def test_driver_uses_pressure_temperature_and_airmass_for_tracer_increment(cell):
    initial = np.full((1, 1, 5), 2.0)
    driver_cell = {k: v for k, v in cell.items() if k != "rho_a"}
    result = fengsha_rrfs_driver(
        **driver_cell,
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[2.0]]),
        dt=10.0,
        dust=initial,
    )
    flux = dust_emission(**{**cell, "rho_a": np.array([[1.2]])})
    np.testing.assert_allclose(result.rate, flux)
    np.testing.assert_allclose(result.bems, flux * 1e9)
    airmas = 1000.0 * 2.0 / GRAV
    pre = initial * 1e-9 + flux * 2.0 * 10.0 / airmas
    expected = pre.copy()
    expected[..., 0] += 0.286 * pre[..., 1]
    expected[..., 4] = 0.714 * pre[..., 1] + pre[..., 2] + pre[..., 3]
    np.testing.assert_allclose(result.dust, expected * 1e9)
    np.testing.assert_allclose(result.emis_dust, result.bems[..., :4].sum(axis=-1))


def test_driver_skipped_cell_keeps_nonzero_initial_tracers(cell):
    initial = np.arange(1.0, 6.0).reshape(1, 1, 5)
    driver_cell = {k: v for k, v in cell.items() if k != "rho_a"}
    result = fengsha_rrfs_driver(
        **{**driver_cell, "xland": np.array([[2.0]])},
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[1.0]]),
        dt=10.0,
        dust=initial,
    )
    np.testing.assert_array_equal(result.dust, initial)
    assert not result.bems.any()
    assert not result.emis_dust.any()


def test_driver_snow_gate_preserves_initial_tracers_without_remap(cell):
    initial = np.arange(1.0, 6.0).reshape(1, 1, 5)
    driver_cell = {k: v for k, v in cell.items() if k != "rho_a"}
    result = fengsha_rrfs_driver(
        **{**driver_cell, "snowh": np.array([[0.001]])},
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[2.0]]),
        dt=10.0,
        dust=initial,
    )
    np.testing.assert_array_equal(result.dust, initial)
    assert not result.rate.any()
    assert not result.bems.any()


def test_driver_zero_drag_cell_skips_tracer_remap(cell):
    """rdrag <= 0 cycles the cell (F:174), so F:369-370 must not remap tracers."""
    initial = np.arange(1.0, 6.0).reshape(1, 1, 5)
    driver_cell = {k: v for k, v in cell.items() if k != "rho_a"}
    result = fengsha_rrfs_driver(
        **{**driver_cell, "rdrag": np.array([[0.0]])},
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[2.0]]),
        dt=10.0,
        dust=initial,
    )
    np.testing.assert_array_equal(result.dust, initial)


def test_driver_tracer_increment_scales_with_dt(cell):
    """The mass added to the tracers is emit * dxy * dt / airmas (F:359,363)."""
    initial = np.zeros((1, 1, 5))
    driver_cell = {k: v for k, v in cell.items() if k != "rho_a"}
    common = dict(
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[2.0]]),
        dust=initial,
    )
    r10 = fengsha_rrfs_driver(**driver_cell, dt=10.0, **common)
    r30 = fengsha_rrfs_driver(**driver_cell, dt=30.0, **common)
    assert r10.dust.any()
    np.testing.assert_allclose(r30.dust, 3.0 * r10.dust, rtol=1e-12)
    np.testing.assert_allclose(r30.bems, r10.bems, rtol=1e-12)


def test_monthly_driver_uses_the_given_month_of_dust12m(cell):
    """The wrapper reads dust12m at jdate(2) (W:298, W:818-822), not a fixed month."""
    driver_cell = {
        k: v for k, v in cell.items() if k not in ("rho_a", "clay", "rdrag", "sand", "ssm", "uthrs")
    }
    host = dict(
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        area=np.array([[2.0]]),
        dt=10.0,
        dust=np.zeros((1, 1, 5)),
    )
    dust12m = np.empty((1, 1, 12, 5))
    for m in range(12):
        dust12m[0, 0, m] = [0.10 + 0.01 * m, 0.90, 0.50, 0.80, 0.10 + 0.005 * m]
    for month in (1, 7, 12):
        fields = {
            name: dust12m[..., month - 1, k]
            for k, name in enumerate(("clay", "rdrag", "sand", "ssm", "uthrs"))
        }
        expected = fengsha_rrfs_driver(**driver_cell, **fields, **host)
        result = fengsha_rrfs_monthly_driver(dust12m, month, **driver_cell, **host)
        assert expected.rate.any()
        np.testing.assert_array_equal(result.rate, expected.rate)
    jan = fengsha_rrfs_monthly_driver(dust12m, 1, **driver_cell, **host)
    jul = fengsha_rrfs_monthly_driver(dust12m, 7, **driver_cell, **host)
    assert not np.allclose(jan.rate, jul.rate)


@pytest.mark.parametrize("month", [0, 13])
def test_monthly_driver_rejects_month_outside_1_to_12(cell, month):
    """A month outside 1-12 reads past the 12-month array in the Fortran."""
    with pytest.raises(ValueError):
        select_month(np.zeros((1, 1, 12, 5)), month)


# --- scaling laws ---


def test_flux_scales_linearly_with_alpha(cell):
    F1 = dust_emission(**cell)
    F2 = dust_emission(
        **{**cell, "settings": replace(cell["settings"], alpha=2.0 * cell["settings"].alpha)}
    )
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-12)


def test_flux_scales_linearly_with_rho_a(cell):
    F1 = dust_emission(**cell)
    F2 = dust_emission(**{**cell, "rho_a": cell["rho_a"] * 3.0})
    assert np.isclose(F2.sum(), F1.sum() * 3.0, rtol=1e-12)


def test_moister_soil_reduces_flux(cell):
    """Higher soil moisture raises the threshold, lowering flux."""
    F_dry = dust_emission(**{**cell, "slc": np.array([[0.01]])})
    F_moist = dust_emission(**{**cell, "slc": np.array([[0.12]])})
    assert F_moist.sum() < F_dry.sum()


# --- source-verified parameter guards ---


def test_kvhmax_value():
    """kvhmax must be 2.0e-4 (source_dust, L307)."""
    assert KVHMAX_DEFAULT == 2.0e-4


def test_rhosoil_is_dead_argument():
    """
    RHOSOIL (2650) is passed to DustEmissionFENGSHA as rhop but never used in
    its body; RHO_DRY (1700) is the live density (Fecan vol->grav). This guard
    documents that the two source constants differ and must not be conflated.
    """
    assert RHOSOIL == 2650.0
    assert RHO_DRY == 1700.0
    assert RHOSOIL != RHO_DRY


def test_mb95_ratio_clay_branches():
    """clay>0.2 -> kvhmax; clay<=0.2 -> 10^(13.4*clay-6)."""
    assert np.isclose(flux_ratio_mb95(np.array(0.3)), KVHMAX_DEFAULT)
    assert np.isclose(flux_ratio_mb95(np.array(0.1)), 10.0 ** (13.4 * 0.1 - 6.0))


def test_uthres_table_values():
    """Spot-check the dry threshold table (dust_data_mod.F90, L51-L66)."""
    assert UTHRES_TABLE[0] == 0.065  # Sand
    assert UTHRES_TABLE[10] == 0.50  # Silty Clay
    assert UTHRES_TABLE[12] == 9999.0  # Other / nonerodible


def test_kok_fractions_sum_to_one():
    assert np.isclose(KOK_FRAC.sum(), 1.0, rtol=1e-12)
    assert (KOK_FRAC > 0.0).all()


# --- initialize() ---


def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("missing-isltyp.npy")


def test_initialize_valid(tmp_path):
    shp = (5, 4)
    p = tmp_path / "isltyp.npy"
    np.save(p, np.full(shp, 3))
    assert initialize(p).shape == shp


# ---------------------------------------------------------------------------
# Host tracer steps of the wrapper (W:796, W:589-590, W:1066-1067)
# ---------------------------------------------------------------------------
def _quiet_step(dust, coarse_pm, **changes):
    """Step on a water cell the driver skips (F:99), so only the host steps act."""
    case = {**_host_step_inputs(dust, coarse_pm), "land": np.array([[2]]), **changes}
    return rrfs_dust_step(**case)


def test_step_floors_a_zero_input_tracer_at_epsilc_so_the_bin_never_starts_at_zero(monkeypatch):
    # W:1066-1067: chem = max(epsilc, gq0). A 0 tracer must enter the step as 1e-12.
    seen = {}
    driver = rrfs.fengsha_rrfs_driver

    def capture(**kwargs):
        seen.update(kwargs)
        return driver(**kwargs)

    monkeypatch.setattr(rrfs, "fengsha_rrfs_driver", capture)
    _quiet_step(0.0, 0.0)
    assert seen["dust"][0, 0, 0] == 1.0e-12
    assert seen["dust"][0, 0, 4] == 1.0e-12


def test_step_bins_two_to_four_start_at_zero_and_bin_five_starts_from_coarse_pm(monkeypatch):
    # W:796 zeroes chem; p_dust_5 = p_coarse_pm (config:69). Old 5-bin code carried any values.
    seen = {}
    driver = rrfs.fengsha_rrfs_driver

    def capture(**kwargs):
        seen.update(kwargs)
        return driver(**kwargs)

    monkeypatch.setattr(rrfs, "fengsha_rrfs_driver", capture)
    rrfs_dust_step(**_host_step_inputs(4.0, 7.0))
    np.testing.assert_array_equal(seen["dust"][0, 0], [4.0, 0.0, 0.0, 0.0, 7.0])


def test_step_emission_does_not_depend_on_old_values_in_bins_two_to_four():
    # Direct check against a 5-bin run: the step equals the driver started from [d,0,0,0,c].
    case = _host_step_inputs(4.0, 7.0)
    step = rrfs_dust_step(**case)
    assert step.rate.sum() > 0
    soil = rrfs.select_month(case["dust12m"], 1)
    ref = fengsha_rrfs_driver(
        **soil,
        ustar=case["ustar"],
        slc=np.array([[0.02]]),
        stemp=np.array([[290.0]]),
        znt=np.array([[0.01]]),
        snowh=np.array([[0.0]]),
        xland=np.array([[1.0]]),
        p8w1=np.array([[100000.0]]),
        p8w2=np.array([[99000.0]]),
        p=np.array([[103338.0]]),
        temperature=np.array([[300.0]]),
        area=case["garea"],
        isltyp=case["soiltyp"],
        dt=10.0,
        dust=np.array([[[4.0, 0.0, 0.0, 0.0, 7.0]]]),
    )
    np.testing.assert_allclose(step.dust, ref.dust[..., 0])
    np.testing.assert_allclose(step.coarse_pm, ref.dust[..., 4])


def test_step_coarse_pm_tracer_is_bin_five_and_survives_a_skipped_cell():
    # Bin 5 is the coarse-PM slot: in a skipped cell it comes back unchanged.
    out = _quiet_step(4.0, 7.0)
    assert out.dust.item() == 4.0
    assert out.coarse_pm.item() == 7.0


def test_step_active_cell_replaces_coarse_pm_by_the_remap_as_the_source_does():
    # F:370 assigns tc(5) = 0.714*tc(2)+tc(3)+tc(4): the old bin 5 is not added back.
    out = rrfs_dust_step(**_host_step_inputs(4.0, 7.0))
    assert out.coarse_pm.item() != 7.0


def test_step_caps_output_tracers_at_200_and_5000():
    # W:589-590: gq0 = min(cap, ...). Too-large tracers are cut, not passed to the host.
    out = _quiet_step(300.0, 6000.0)
    assert out.dust.item() == 200.0
    assert out.coarse_pm.item() == 5000.0
    below = _quiet_step(199.0, 4999.0)
    assert (below.dust.item(), below.coarse_pm.item()) == (199.0, 4999.0)


def test_step_floors_output_tracers_at_epsilc():
    # W:589-590: max(epsilc, ...). In an active cell with no emission F:370 sets bin 5 to 0;
    # the host must still get 1e-12, never 0.
    out = rrfs_dust_step(**_host_step_inputs(4.0, 7.0, ustar=0.0))
    assert out.rate.sum() == 0
    assert out.coarse_pm.item() == 1.0e-12


def test_step_setting_names_carry_a_source_line():
    meta = {f.name: f.metadata for f in rrfs.fields(RrfsSettings)}
    assert meta["epsilc"]["source"] == "rrfs_smoke_config.F90:17"
    assert meta["dust_cap"]["source"] == "rrfs_smoke_wrapper.F90:589"
    assert meta["coarse_pm_cap"]["source"] == "rrfs_smoke_wrapper.F90:590"
    assert (RrfsSettings().epsilc, RrfsSettings().dust_cap, RrfsSettings().coarse_pm_cap) == (
        1.0e-12,
        200.0,
        5000.0,
    )
