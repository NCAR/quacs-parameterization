"""Tests for quacs.dust.fengsha_ls99."""

import re
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import quacs.dust.fengsha_ls99.fengsha_ls99 as ls99
from quacs.dust.fengsha_ls99.fengsha_ls99 import (
    A_THRESH,
    B_THRESH,
    CFAC,
    DP,
    DUST_DENS_CMAQ,
    FRACMJ,
    FRACMK,
    GPKG,
    GRAV,
    HB_LAMBDAB,
    HV_SEASONAL,
    LAMBDAB,
    LANDUSE_SCHEMES,
    PI,
    SOILML1,
    SOILTXT,
    TINY_REAL4,
    DustInputWarning,
    Ls99Settings,
    _get_hv,
    _host_landuse,
    _remap_soil_type,
    _surface_rate,
    dust_emission,
    dust_hflux,
    fengsha_dust_driver,
    tunable_names,
)


# ============================================================================
# Embedded pytest suite
# ============================================================================
def test_soiltxt_matches_source():
    """SOILTXT rows must match DUST_EMIS.F's soiltxt literal exactly."""
    expected = {
        4: [0.00, 0.17, 0.70, 0.13],
        5: [0.00, 0.10, 0.85, 0.05],
        6: [0.00, 0.43, 0.39, 0.18],
        7: [0.29, 0.29, 0.15, 0.27],
        8: [0.00, 0.10, 0.56, 0.34],
        9: [0.00, 0.32, 0.34, 0.34],
        10: [0.00, 0.52, 0.06, 0.42],
        11: [0.00, 0.06, 0.47, 0.47],
    }
    for idx, vals in expected.items():
        assert np.allclose(SOILTXT[idx - 1, :], vals), f"SOILTXT row {idx} mismatch"


def test_soil_type_remap():
    """Remap must match DUST_EMIS.F lines 829-833."""
    assert _remap_soil_type(np.array(1)) == 1
    assert _remap_soil_type(np.array(4)) == 4
    assert _remap_soil_type(np.array(5)) == 6  # j>4 -> j+1
    assert _remap_soil_type(np.array(12)) == 13  # j>4 -> j+1
    assert _remap_soil_type(np.array(13)) == 13  # j>4 -> 14, then >13 -> 13
    assert _remap_soil_type(np.array(16)) == 13


def test_dust_hflux_threshold_gating():
    """Below threshold u*, hflux must be exactly zero."""
    hflux_low = dust_hflux(DP, SOILTXT[0], 1.0, 1.0, ustr=0.01, sd_ep=0.1, dens=1.2)
    assert hflux_low == 0.0


def test_dust_hflux_positive_above_threshold():
    hflux_high = dust_hflux(DP, SOILTXT[0], 1.0, 1.0, ustr=5.0, sd_ep=0.1, dens=1.2)
    assert hflux_high > 0.0


def test_landuse_scheme_registry_shapes():
    for name, reg in LANDUSE_SCHEMES.items():
        assert len(reg["categories"]) == reg["n_dlcat"]
        assert len(reg["dmap"]) == reg["n_dlcat"]


def test_dust_emission_output_shape_and_nonnegative():
    u10m = np.array([[15.0]])
    soim = np.array([[0.02]])
    dens = np.array([[1.2]])
    snow = np.array([[0.0]])
    rain = np.array([[0.0]])
    veg = np.array([[0.05]])
    stype = np.array([[1]])
    luf = np.zeros((1, 1, 4))
    luf[0, 0, 1] = 80.0  # USGS24 col index 1 = barren

    F = dust_emission(u10m, soim, dens, snow, rain, veg, stype, luf, jday=100, scheme="usgs24")
    assert F.shape == (1, 1, 2)
    assert (F >= 0).all()


def test_dust_emission_snow_suppresses():
    u10m = np.array([[15.0]])
    soim = np.array([[0.02]])
    dens = np.array([[1.2]])
    snow = np.array([[0.5]])
    rain = np.array([[0.0]])
    veg = np.array([[0.05]])
    stype = np.array([[1]])
    luf = np.zeros((1, 1, 4))
    luf[0, 0, 1] = 80.0

    F = dust_emission(u10m, soim, dens, snow, rain, veg, stype, luf, jday=100, scheme="usgs24")
    assert F[0, 0, :].sum() == 0.0


def test_dust_emission_calm_wind_near_zero():
    u10m = np.array([[3.0]])
    soim = np.array([[0.02]])
    dens = np.array([[1.2]])
    snow = np.array([[0.0]])
    rain = np.array([[0.0]])
    veg = np.array([[0.05]])
    stype = np.array([[1]])
    luf = np.zeros((1, 1, 4))
    luf[0, 0, 1] = 80.0

    F = dust_emission(u10m, soim, dens, snow, rain, veg, stype, luf, jday=100, scheme="usgs24")
    assert F[0, 0, :].sum() < 1.0e-8


def test_mode_split_sums_to_one():
    assert np.isclose(FRACMJ + FRACMK, 1.0)


def test_scheme_consistent_across_schemes_with_equivalent_coverage():
    """Same effective shrub/barren coverage under different schemes should
    give the same flux, since physics is BELD3-class-indexed."""
    u10m = np.array([[15.0]])
    soim = np.array([[0.02]])
    dens = np.array([[1.2]])
    snow = np.array([[0.0]])
    rain = np.array([[0.0]])
    veg = np.array([[0.05]])
    stype = np.array([[1]])

    luf_usgs24 = np.zeros((1, 1, 4))
    luf_usgs24[0, 0, 0] = 80.0  # 80% shrubland (BELD3 class 1)

    luf_nlcd40 = np.zeros((1, 1, 6))
    luf_nlcd40[0, 0, 0] = 80.0  # 80% closed shrublands (also BELD3 class 1)

    F1 = dust_emission(
        u10m, soim, dens, snow, rain, veg, stype, luf_usgs24, jday=100, scheme="usgs24"
    )
    F2 = dust_emission(
        u10m, soim, dens, snow, rain, veg, stype, luf_nlcd40, jday=100, scheme="nlcd40"
    )
    assert np.allclose(F1, F2)


def test_invalid_scheme_raises():
    u10m = np.array([[15.0]])
    soim = np.array([[0.02]])
    dens = np.array([[1.2]])
    snow = np.array([[0.0]])
    rain = np.array([[0.0]])
    veg = np.array([[0.05]])
    stype = np.array([[1]])
    luf = np.zeros((1, 1, 4))
    try:
        dust_emission(
            u10m, soim, dens, snow, rain, veg, stype, luf, jday=100, scheme="not_a_real_scheme"
        )
        assert False, "expected ValueError"
    except ValueError:
        pass


def _port_case():
    def a(value):
        return np.array([[value]], dtype=float)

    luf = np.zeros((1, 1, 24))
    luf[0, 0, 18] = 1.0
    return dict(
        u10m=a(30),
        soim_vol=a(0),
        dens=np.array([[[1.1], [2.2]]]),
        snocov=a(0),
        rn=a(0),
        rc=a(0),
        vegfrac=a(0.1),
        soil_type=a(2),
        luf=luf,
        jdate=2000140,
        lwmask=a(1),
        layer_height=a(50),
        jacobm=a(1),
        wmax_last=a(0),
        sd_ep_last=a(0),
        settings=Ls99Settings(spcfac=((1.0, 0.0), (0.5, 0.25), (0.0, 1.0))),
    )


def test_explicit_default_settings_equal_implicit_defaults():
    x = _port_case()
    direct = dict(
        u10m=x["u10m"],
        soim_vol=x["soim_vol"],
        dens=x["dens"][:, 0, :],
        snocov=x["snocov"],
        rn=x["rn"],
        rc=x["rc"],
        vegfrac=x["vegfrac"],
        soil_type=x["soil_type"],
        luf=np.array([[[0.0, 100.0, 0.0, 0.0]]]),
        jday=140,
    )
    np.testing.assert_array_equal(
        dust_emission(**direct), dust_emission(**direct, settings=Ls99Settings())
    )


def _scale(value):
    return value * 1.2


def _scale_first(value):
    return (value[0] * 1.2, *value[1:])


def _scale_third(value):
    return (*value[:2], value[2] * 1.2, value[3])


def _scale_one_month(value):
    changed = list(value)
    month = list(changed[3])
    month[2] *= 1.2
    changed[3] = tuple(month)
    return tuple(changed)


def _scale_one_texture(value):
    changed = list(value)
    changed[1] = (value[1][0] * 1.2, *value[1][1:])
    return tuple(changed)


# One entry per setting. A setting without an entry fails the test below.
LS99_SETTING_CHANGES = {
    "rhop": _scale,
    "a_thresh": _scale,
    "b_thresh": _scale,
    "amen": _scale,
    "cfac": _scale,
    "dp": _scale_first,
    "eropot": lambda value: (value[0], value[1] * 1.2, value[2]),
    "soiltxt": _scale_one_texture,
    "soilml1": lambda value: (value[0], 0.01, *value[2:]),
    "flxfac1": _scale_first,
    "flxfac2": _scale_first,
    "kvh_offset": _scale,
    "wmax_clay": _scale,
    "wmax_base": _scale,
    "soimt_base": _scale,
    "soimt_sand": _scale,
    "moist_cutoff": lambda value: 0.1,
    "moist_scale": _scale,
    "moist_power": _scale,
    "mv": _scale,
    "sigv": _scale,
    "betav": _scale,
    "mb": _scale,
    "sigb": _scale,
    "betab": _scale,
    "lambda_v": _scale,
    "veg_max": lambda value: 0.8,
    "veg_min": lambda value: 0.1,
    "lambdab": _scale_third,
    "hb_lambdab": _scale_third,
    "hv_seasonal": _scale_one_month,
    "z0_low_scale": _scale,
    "z0_low_power": _scale,
    "z0_split": lambda value: 0.01,
    "z0_high_scale": _scale,
    "z0_high_power": _scale,
    "fruf_cap": _scale,
    "wind_height": _scale,
    "rain_max": lambda value: -0.001,
    "snow_max": lambda value: -0.001,
    "fracmj": _scale,
    "fracmk": _scale,
    "dust_dens": _scale_first,
    "spcfac": lambda value: ((0.8, 0.0), *value[1:]),
    "dgv": _scale_first,
}
# fruf_cap only acts when the drag factor is capped, so that case needs a large sigb.
LS99_SETTING_SEED = {"fruf_cap": {"sigb": 10000.0}}


def _ls99_outputs(settings):
    original = _port_case()
    original["soim_vol"][:] = 0.08
    cases = [
        original,
        {**original, "vegfrac": np.array([[0.0]])},
        {**original, "vegfrac": np.array([[0.97]])},
        {**original, "soil_type": np.array([[3.0]])},
    ]
    pieces = []
    for case in cases:
        result = fengsha_dust_driver(**{**case, "settings": settings})
        pieces.extend(
            result[name].ravel()
            for name in (
                "dust_em",
                "fmoit",
                "soimt",
                "wmax",
                "fruf",
                "kvh",
                "ustr",
                "mode_mass",
                "number_rates",
                "surface_area_rates",
                "species_rates",
            )
        )
    return np.concatenate(pieces)


@pytest.mark.parametrize("name", tunable_names(Ls99Settings))
def test_each_setting_reaches_an_output_so_tuning_it_is_never_silently_ignored(name):
    seed = replace(_port_case()["settings"], **LS99_SETTING_SEED.get(name, {}))
    changed = LS99_SETTING_CHANGES[name](getattr(seed, name))
    assert not np.array_equal(
        _ls99_outputs(replace(seed, **{name: changed})), _ls99_outputs(seed), equal_nan=True
    ), name


def test_every_setting_case_names_a_real_setting():
    assert set(LS99_SETTING_CHANGES) <= set(tunable_names(Ls99Settings))


def test_each_mode_density_changes_only_its_own_mode():
    """Mode 1 and mode 2 densities feed separate number and area rates, so both defaults matter."""
    case = _port_case()
    settings = case["settings"]
    base = fengsha_dust_driver(**case)
    for mode in (0, 1):
        dens = list(settings.dust_dens)
        dens[mode] *= 1.2
        out = fengsha_dust_driver(**{**case, "settings": replace(settings, dust_dens=tuple(dens))})
        other = 1 - mode
        assert out["number_rates"][..., mode] != base["number_rates"][..., mode]
        assert out["surface_area_rates"][..., mode] != base["surface_area_rates"][..., mode]
        assert out["number_rates"][..., other] == base["number_rates"][..., other]


def test_default_mode_densities_are_the_first_two_cmaq_values():
    assert tuple(float(v) for v in Ls99Settings().dust_dens) == (2200.0, 2156.55)


def test_wind_height_sets_the_roughness_log_in_friction_velocity():
    """ustr = karman * u10 / log(height / z0); a tuned height must reach it."""
    case = _port_case()
    base = fengsha_dust_driver(**case)
    tall = fengsha_dust_driver(**{**case, "settings": replace(case["settings"], wind_height=20.0)})
    assert tall["ustr"][..., 1] != base["ustr"][..., 1]


def test_default_numbers_are_pinned_so_a_default_drift_fails():
    """Defaults are the CMAQ values; a 1e-6 drift in any of them must show here."""
    case = _port_case()
    out = fengsha_dust_driver(**case)
    expected = {
        "dust_em": [0.00030192938627473164],
        "mode_mass": [2.1135057039231218e-05, 0.00028079432923550047],
        "number_rates": [59181714.41701107, 14855322.48742619],
        "surface_area_rates": [5.2675545423201036e-05, 0.00018888867794445964],
        "species_rates": [
            2.1135057039231218e-05,
            1.0567528519615609e-05,
            0.0,
            0.0,
            7.019858230887512e-05,
            0.00028079432923550047,
        ],
        "ustr": [0.0, 1.467627240154721, 0.0, 0.0, 0.0],
        "kvh": [0.0, 0.00014556419116459603, 0.0, 0.0, 0.0],
    }
    assert out["dust_em"].min() > 0
    for name, values in expected.items():
        np.testing.assert_allclose(
            np.asarray(out[name]).ravel(), values, rtol=1e-9, atol=0, err_msg=name
        )
    case["soim_vol"][:] = 0.08
    wet = fengsha_dust_driver(**case)
    np.testing.assert_allclose(wet["dust_em"].ravel(), [0.00020099780818433582], rtol=1e-9)
    np.testing.assert_allclose(
        wet["number_rates"].ravel(), [39397936.81290324, 9889356.238836654], rtol=1e-9
    )


def test_density_and_species_fractions_leave_mass_flux_unchanged():
    case = _port_case()
    base = fengsha_dust_driver(**case)
    settings = case["settings"]
    dense = fengsha_dust_driver(
        **{
            **case,
            "settings": replace(
                settings, dust_dens=(settings.dust_dens[0] * 1.2, settings.dust_dens[1])
            ),
        }
    )
    fractions = fengsha_dust_driver(
        **{**case, "settings": replace(settings, spcfac=((0.8, 0.0), *settings.spcfac[1:]))}
    )
    np.testing.assert_array_equal(dense["dust_em"], base["dust_em"])
    np.testing.assert_array_equal(fractions["dust_em"], base["dust_em"])
    assert not np.array_equal(dense["number_rates"], base["number_rates"])
    assert not np.array_equal(dense["surface_area_rates"], base["surface_area_rates"])
    assert not np.array_equal(fractions["species_rates"], base["species_rates"])
    np.testing.assert_array_equal(fractions["number_rates"], base["number_rates"])
    np.testing.assert_array_equal(fractions["surface_area_rates"], base["surface_area_rates"])


def test_cmaq_gravity_and_rounded_threshold_literals_keep_flux_pinned():
    assert GRAV == 9.80622
    assert A_THRESH == 260.60061
    assert B_THRESH == 1.6540342e-6
    assert CFAC == 1000.0 / 9.80622


def test_land_water_mask_prevents_water_cell_emissions():
    x = _port_case()
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] > 0
    x["lwmask"][:] = 0
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0


def test_march_begins_after_day_59():
    assert _get_hv(59) is HV_SEASONAL["nov_feb"]
    assert _get_hv(60) is HV_SEASONAL["mar"]
    assert _get_hv(61) is HV_SEASONAL["mar"]


def test_reconstructed_silt_fraction_reaches_horizontal_flux(monkeypatch):
    seen = []

    def capture(dp, texture, *args):
        seen.append(texture.copy())
        return 0.0

    monkeypatch.setattr(ls99, "dust_hflux", capture)
    x = _port_case()
    x["soil_type"][:] = 13
    x["soil_type_is_premapped"] = True
    fengsha_dust_driver(**x)
    assert seen and np.array_equal(seen[0], [0.0, 0.0, 1.0, 0.0])


def test_px_texture_replaces_lookup_and_all_16_rows_exist(monkeypatch):
    assert len(SOILTXT) == len(SOILML1) == 16
    assert SOILML1[15] == 0.280
    assert np.array_equal(SOILTXT[11], [0.0, 0.22, 0.20, 0.58])
    seen = []

    def capture(dp, texture, *args):
        seen.append(texture.copy())
        return 0.0

    monkeypatch.setattr(ls99, "dust_hflux", capture)
    x = _port_case()
    x.update(csand=np.array([[0.05]]), fmsand=np.array([[0.15]]), clay=np.array([[0.55]]))
    fengsha_dust_driver(**x)
    assert seen and np.allclose(seen[0], [0.05, 0.15, 0.25, 0.55])
    x["soil_type"][:] = 16
    x["soil_type_is_premapped"] = True
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0


def test_winter_barren_ustar_uses_computed_z0_below_one_ten_thousand():
    x = _port_case()
    x["vegfrac"][:] = 0
    x["jdate"] = 2000001
    out = fengsha_dust_driver(**x)
    vf = 0.005
    lv = -0.35 * np.log(1 - vf)
    lam = LAMBDAB[3] + lv
    height = (HB_LAMBDAB[3] + _get_hv(1)[3] * lv) / lam
    z0 = 0.96 * lam**1.07 * height
    assert z0 < 1e-4
    assert np.isclose(out["ustr"][0, 0, 1], 0.4 * 30 / np.log(10 / z0), rtol=1e-14, atol=0)


def test_bad_luf_length_and_out_of_range_elus_stop():
    x = _port_case()
    x["luf"] = np.zeros((1, 1, 3))
    with pytest.raises(ValueError, match="luf"):
        fengsha_dust_driver(**x)
    x = _port_case()
    x["luf"][0, 0, 18] = 1.01 / 0.9
    with pytest.raises(ValueError, match="elus"):
        fengsha_dust_driver(**x)


def test_qam_resets_on_each_call_after_an_emitting_cell():
    x = _port_case()
    assert fengsha_dust_driver(**x)["qam"].sum() > 0
    x["luf"][:] = 0
    assert np.array_equal(fengsha_dust_driver(**x)["qam"], np.zeros((1, 1, 5)))


def test_impossible_nonfatal_input_uses_one_warning_class():
    x = _port_case()
    x["u10m"][:] = -1
    with pytest.warns(DustInputWarning):
        out = fengsha_dust_driver(**x)
    assert np.isfinite(out["dust_em"]).all()


def test_two_rain_fields_and_snow_limit_gate_emissions():
    x = _port_case()
    for name in ("rn", "rc"):
        x[name][:] = 0.026
        assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0
        x[name][:] = 0
    x["rn"][:] = 0.0127
    x["rc"][:] = 0.0127
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] > 0
    x["snocov"][:] = 0.001
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0


def test_driver_qam_species_moments_and_diagnostic_totals_share_rate():
    x = _port_case()
    fac = np.asarray(x["settings"].spcfac)
    out = fengsha_dust_driver(**x, diagnostic_totals=5.0, sync_seconds=30.0)
    assert np.isclose(out["qam"].sum(), out["dust_em"].sum())
    assert np.isclose(out["dust_em"][0, 0] * 50, out["surface_rate"][0, 0])
    assert np.allclose(out["species_rates"][0, 0], out["mode_mass"][0, 0, :, None] * fac.T)
    assert (out["number_rates"] > 0).all() and (out["surface_area_rates"] > 0).all()
    assert np.allclose(out["diagnostic_totals"], 5 + 30 * out["diagnostics"])
    assert np.isclose(out["diagnostics"][0, 0, 5], out["dust_em"][0, 0])


def test_saved_wmax_and_sd_ep_survive_a_gated_call():
    x = _port_case()
    first = fengsha_dust_driver(**x)
    assert first["wmax"][0, 0] >= 0 and first["sd_ep"][0, 0] > 0
    x["wmax_last"] = first["wmax"]
    x["sd_ep_last"] = first["sd_ep"]
    x["lwmask"][:] = 0
    second = fengsha_dust_driver(**x)
    assert np.array_equal(second["wmax"], first["wmax"])
    assert np.array_equal(second["sd_ep"], first["sd_ep"])
    assert second["dust_em"][0, 0] == 0


def test_rjacm_is_the_reciprocal_of_jacobm():
    # CMAQ reads JACOBM and sets RJACM = 1 / JACOBM (ASX_DATA_MOD.F:878-882).
    # The layer rate qam carries rjacm / layer_height, the surface rate does not.
    x = _port_case()
    one = fengsha_dust_driver(**x)
    x["jacobm"][:] = 0.5
    half = fengsha_dust_driver(**x)
    assert np.allclose(half["qam"], 2 * one["qam"], rtol=1e-12)
    assert np.allclose(half["surface_rate"], one["surface_rate"], rtol=1e-12)


def test_mpas_sets_rjacm_to_one_whatever_jacobm_is():
    # ASX_DATA_MOD.F:875-876 under #ifdef mpas: RJACM = 1.0, JACOBM is never read.
    x = _port_case()
    plain = fengsha_dust_driver(**x)
    x["soim_vol"] = np.zeros((1, 2, 1))
    reference = fengsha_dust_driver(**x, mpas=True)
    assert np.array_equal(reference["qam"], plain["qam"])
    x["jacobm"][:] = 0.25
    assert np.array_equal(fengsha_dust_driver(**x, mpas=True)["qam"], reference["qam"])
    x["jacobm"] = None
    assert np.array_equal(fengsha_dust_driver(**x, mpas=True)["qam"], reference["qam"])


def test_jacobm_is_required_without_mpas():
    # Without mpas CMAQ must read JACOBM; silently using 1.0 would hide a missing field.
    x = _port_case()
    x["jacobm"] = None
    with pytest.raises(ValueError, match="jacobm"):
        fengsha_dust_driver(**x)


def test_default_start_values_equal_zero_fields():
    # CMAQ sets wmax and sd_ep to 0 at start (DUST_EMIS.F:281-282). A gated cell
    # returns its start values unchanged, so it shows what the default was.
    x = _port_case()
    x["lwmask"][:] = 0
    default = _port_case()
    default["lwmask"][:] = 0
    del default["wmax_last"], default["sd_ep_last"]
    out = fengsha_dust_driver(**default)
    zeros = fengsha_dust_driver(**x)
    assert np.array_equal(out["wmax"], np.zeros((1, 1)))
    assert np.array_equal(out["sd_ep"], np.zeros((1, 1)))
    assert np.array_equal(out["diagnostics"], zeros["diagnostics"])
    y = _port_case()
    del y["wmax_last"], y["sd_ep_last"]
    assert np.array_equal(
        fengsha_dust_driver(**y)["qam"], fengsha_dust_driver(**_port_case())["qam"]
    )


def test_later_step_with_saved_fields_differs_from_a_fresh_start_only_where_cmaq_keeps_state():
    x = _port_case()
    first = fengsha_dust_driver(**x)
    x["wmax_last"], x["sd_ep_last"] = first["wmax"], first["sd_ep"]
    fresh = {k: v for k, v in x.items() if k not in ("wmax_last", "sd_ep_last")}
    # Emitting cell: both fields are recomputed, so the saved ones change nothing.
    assert np.array_equal(fengsha_dust_driver(**x)["qam"], fengsha_dust_driver(**fresh)["qam"])
    # Gated cell: CMAQ keeps the saved fields; a fresh start loses them.
    x["lwmask"][:] = 0
    fresh["lwmask"][:] = 0
    kept, lost = fengsha_dust_driver(**x), fengsha_dust_driver(**fresh)
    assert np.array_equal(kept["sd_ep"], first["sd_ep"]) and kept["sd_ep"][0, 0] > 0
    assert lost["sd_ep"][0, 0] == 0
    assert not np.array_equal(kept["diagnostics"], lost["diagnostics"])


def test_zero_mass_moments_use_cmaq_real4_tiny_and_pinned_pi():
    x = _port_case()
    x["luf"][:] = 0
    out = fengsha_dust_driver(**x)
    factnum = 1e18 * np.exp(4.5 * np.log(2.0) ** 2) / np.array([1.3914, 5.2590]) ** 3
    assert PI == 3.14159265 and GPKG == 1e3
    assert np.array_equal(DUST_DENS_CMAQ, [2200.0, 2156.55, 2536.92])
    assert np.allclose(out["number_rates"][0, 0], TINY_REAL4 * factnum, rtol=1e-15, atol=0)


def test_missing_species_or_densities_raise_before_flux_calculation():
    x = _port_case()
    x["settings"] = replace(x["settings"], spcfac=None)
    with pytest.raises(ValueError, match="spcfac and dust_dens"):
        fengsha_dust_driver(**x)
    x = _port_case()
    x["settings"] = replace(x["settings"], dust_dens=None)
    with pytest.raises(ValueError, match="spcfac and dust_dens"):
        fengsha_dust_driver(**x)


def test_partial_px_texture_is_rejected():
    for name in ("csand", "fmsand"):
        x = _port_case()
        x[name] = np.array([[0.2]])
        with pytest.raises(ValueError, match="PX texture"):
            fengsha_dust_driver(**x)


def test_modis_placeholder_has_no_desert_mapping():
    x = _port_case()
    x["scheme"] = "modis"
    x["luf"] = np.zeros((1, 1, 1))
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0
    x["luf"][0, 0, 0] = 1.0
    assert fengsha_dust_driver(**x)["dust_em"][0, 0] == 0


def test_driver_extracts_jday_across_month_boundary():
    x = _port_case()
    x["jdate"] = 2023059
    winter = fengsha_dust_driver(**x)["dust_em"]
    x["jdate"] = 2023060
    march = fengsha_dust_driver(**x)["dust_em"]
    assert winter.sum() > 0 and march.sum() > 0
    assert not np.array_equal(winter, march)


def test_driver_uses_first_density_level():
    x = _port_case()
    first = fengsha_dust_driver(**x)["dust_em"]
    x["dens"][:, 1, :] = 8.0
    np.testing.assert_array_equal(fengsha_dust_driver(**x)["dust_em"], first)
    x["dens"][:, 0, :] = 1.5
    assert not np.array_equal(fengsha_dust_driver(**x)["dust_em"], first)


def test_driver_mpas_first_soil_layer_and_non_mpas_ready_value():
    x = _port_case()
    dry = fengsha_dust_driver(**x)["dust_em"]
    x["soim_vol"] = np.array([[[0.0], [0.9]]])
    np.testing.assert_array_equal(fengsha_dust_driver(**x, mpas=True)["dust_em"], dry)
    x["soim_vol"][:, 0, :] = 0.5
    assert fengsha_dust_driver(**x, mpas=True)["dust_em"].sum() < dry.sum()


def test_driver_clamps_negative_snow_with_warning(monkeypatch):
    x = _port_case()
    x["snocov"][:] = -0.5
    original = _surface_rate
    seen = []

    def capture(*args):
        seen.append(np.array(args[3], copy=True))
        return original(*args)

    monkeypatch.setattr(ls99, "_surface_rate", capture)
    with pytest.warns(DustInputWarning, match="snocov"):
        negative = fengsha_dust_driver(**x)["dust_em"]
    np.testing.assert_array_equal(seen[0], 0.0)
    x["snocov"][:] = 0
    np.testing.assert_array_equal(negative, fengsha_dust_driver(**x)["dust_em"])


def test_driver_selects_desert_column_and_converts_fraction_to_percent():
    x = _port_case()
    out = fengsha_dust_driver(**x)
    assert np.isclose(out["elus"][0, 0, 1], 90.0)
    x["luf"][0, 0, 18] = 0
    x["luf"][0, 0, 1] = 1
    assert fengsha_dust_driver(**x)["dust_em"].sum() == 0


def test_driver_desert_maps_for_each_landuse_scheme():
    for scheme, ncat, full_column, desert_column in (
        ("usgs24", 24, 19, 1),
        ("modis_noah", 20, 16, 2),
        ("nlcd_modis", 50, 47, 6),
        ("nlcd40", 40, 27, 2),
    ):
        x = _port_case()
        x["scheme"] = scheme
        x["luf"] = np.zeros((1, 1, ncat))
        x["luf"][0, 0, full_column - 1] = 0.4
        out = fengsha_dust_driver(**x)
        assert np.isclose(out["elus"][0, 0, desert_column], 36.0)
        assert out["dust_em"].sum() > 0


def test_host_desert_columns_match_pinned_lus_map():
    # lus_data_module.F:73-77,158-163,204-212,289-296 (1-based LUFRAC indices).
    expected = {
        "usgs24": (24, (8, 19, 22, 23)),
        "modis_noah": (20, (6, 7, 16, 20)),
        "nlcd_modis": (50, (37, 38, 7, 12, 13, 18, 47)),
        "nlcd40": (40, (6, 7, 27, 31, 32, 16)),
    }
    for scheme, (ncat, columns) in expected.items():
        for column in range(1, ncat + 1):
            full = np.zeros((1, 1, ncat))
            full[0, 0, column - 1] = 0.2
            desert, _ = _host_landuse(full, scheme, (1, 1))
            wanted = np.zeros(len(columns))
            if column in columns:
                wanted[columns.index(column)] = 20.0
            np.testing.assert_array_equal(
                desert[0, 0], wanted, err_msg=f"{scheme} LUFRAC_{column:02d}"
            )


def test_host_urban_classes_match_pinned_lus_groups():
    # LUS_DEFN.F:128,176-179,227-230,270.
    expected = {
        "usgs24": (24, (1,)),
        "modis_noah": (20, (13,)),
        "nlcd_modis": (50, (4, 5, 6, 44)),
        "nlcd40": (40, (24, 25, 26, 13)),
    }
    for scheme, (ncat, columns) in expected.items():
        for column in range(1, ncat + 1):
            full = np.zeros((1, 1, ncat))
            full[0, 0, column - 1] = 0.2
            _, cover = _host_landuse(full, scheme, (1, 1))
            assert cover[0, 0, 2] == (20.0 if column in columns else 0.0), (scheme, column)


def test_host_forest_classes_match_pinned_lus_groups():
    # LUS_DEFN.F:131-138,182-194,233-245,273-278.
    expected = {
        "usgs24": (24, (11, 12, 13, 14, 15, 18, 21)),
        "modis_noah": (20, (1, 2, 3, 4, 5, 11)),
        "nlcd_modis": (50, (9, 10, 11, 21, 22, 24, 32, 33, 34, 35, 36, 39, 42)),
        "nlcd40": (40, (28, 29, 30, 39, 40, 1, 2, 3, 4, 5, 8, 11)),
    }
    for scheme, (ncat, columns) in expected.items():
        for column in range(1, ncat + 1):
            full = np.zeros((1, 1, ncat))
            full[0, 0, column - 1] = 0.2
            _, cover = _host_landuse(full, scheme, (1, 1))
            assert cover[0, 0, 3] == (20.0 if column in columns else 0.0), (scheme, column)


def test_host_cover_multiplies_each_fraction_before_fortran_order_sum():
    # LUS_DEFN.F:99,176-179: NLCD40 urban terms are added in this order.
    values = (
        0.020172357786591277,
        0.010141581433968743,
        0.0036367025695007186,
        0.028801469951402436,
    )
    full = np.zeros((1, 1, 40))
    for column, value in zip((24, 25, 26, 13), values):
        full[0, 0, column - 1] = value
    _, cover = _host_landuse(full, "nlcd40", (1, 1))
    expected = sum(value * 100.0 for value in values)
    assert cover[0, 0, 2] == expected


def test_host_cover_rescale_uses_pinned_real4_factor():
    # LUS_DEFN.F:298 casts 100/tland to kind 4 before multiplying each group.
    full = np.zeros((1, 1, 24))
    full[0, 0, [18, 0, 10]] = [1.0, 1.0, 0.75]
    _, cover = _host_landuse(full, "usgs24", (1, 1))
    assert cover[0, 0, 2] == 100.0 * float(np.float32(100.0 / 275.0))


def test_driver_merges_clamps_and_rescales_urban_forest_diagnostics():
    x = _port_case()
    x["luf"][0, 0, 0] = 1.5
    x["luf"][0, 0, 10] = 0.75
    out = fengsha_dust_driver(**x)
    n = LANDUSE_SCHEMES["usgs24"]["n_dlcat"]
    urban = out["diagnostics"][0, 0, 5 * (n + 1) + 6]
    forest = out["diagnostics"][0, 0, 5 * (n + 1) + 7]
    assert np.isclose(urban, 100 * 100 / 275)
    assert np.isclose(forest, 75 * 100 / 275)


def test_driver_merges_multiple_urban_and_forest_classes():
    x = _port_case()
    x["luf"][0, 0, 18] = 0
    x["scheme"] = "nlcd_modis"
    x["luf"] = np.zeros((1, 1, 50))
    x["luf"][0, 0, [3, 4, 8, 9]] = [0.1, 0.2, 0.15, 0.25]
    out = fengsha_dust_driver(**x)
    n = LANDUSE_SCHEMES["nlcd_modis"]["n_dlcat"]
    assert np.isclose(out["diagnostics"][0, 0, 5 * (n + 1) + 6], 30.0)
    assert np.isclose(out["diagnostics"][0, 0, 5 * (n + 1) + 7], 40.0)


AERO_DATA = Path(__file__).resolve().parents[1] / "sources" / "CCTM__src__aero__aero6__AERO_DATA.F"


# Copy of the CMAQ dust_spc table, AERO_DATA.F:505-532, as written in the source:
# (name, molecular weight, fine fraction, coarse fraction). The test below checks this
# copy against the pinned file, so the copy guards the Python list where the file is absent.
SOURCE_DUST_SPC = (
    ("SO4", 96.0, 0.02250, 0.02655),
    ("NO3", 62.0, 0.00020, 0.00160),
    ("CL", 35.5, 0.00945, 0.01190),
    ("NH4", 18.0, 0.00005, 0.0),
    ("NA", 23.0, 0.03935, 0.0),
    ("CA", 40.1, 0.07940, 0.0),
    ("MG", 24.3, 0.01900, 0.0),
    ("K", 39.1, 0.03770, 0.0),
    ("POC", 220.0, 0.01075, 0.0),
    ("PNCOM", 220.0, 0.00430, 0.0),
    ("LVPO1", 218.0, 0.01501, 0.0),
    ("LVOO1", 136.0, 2.23e-5, 0.0),
    ("EC", 12.0, 0.0, 0.0),
    ("FE", 55.8, 0.03355, 0.0),
    ("AL", 27.0, 0.05695, 0.0),
    ("SI", 28.1, 0.19425, 0.0),
    ("TI", 47.9, 0.00280, 0.0),
    ("MN", 54.9, 0.00115, 0.0),
    ("H2O", 18.0, 0.00541, 0.00637),
    ("OTHR", 200.0, 0.48319, 0.0),
    ("SOIL", 100.0, 0.0, 0.95358),
    ("MN_HAPS", 54.9, 0.00041, 0.00041),
    ("NI", 58.7, 2.0e-05, 2.0e-05),
    ("CR_III", 52.0, 5.6e-05, 5.6e-05),
    ("AS", 74.92, 5.5e-06, 5.5e-06),
    ("PB", 207.2, 1.6e-05, 1.6e-05),
    ("CD", 112.4, 9.7e-08, 9.7e-08),
    ("PHG", 200.5, 1.4e-08, 1.4e-08),
)


def test_default_species_table_equals_the_cmaq_copy():
    """The default split must be the CMAQ dust_spc table, row for row.

    The comparison with the compiled Fortran feeds an invented 3-row table to both codes.
    It cannot catch a missing or mistyped default. This test runs without the pinned Fortran files.
    """
    expected = tuple((n, float(mw), float(f), float(c)) for n, mw, f, c in SOURCE_DUST_SPC)
    assert len(expected) == 28
    assert ls99.DUST_SPC == expected
    assert Ls99Settings().spcfac == tuple((f, c) for _, _, f, c in expected)


def test_cmaq_copy_matches_every_pinned_aero_data_row():
    """The copy above must equal the pinned AERO_DATA.F. Runs only where the pin exists."""
    if not AERO_DATA.exists():
        pytest.skip(f"pinned CMAQ file not present: {AERO_DATA}")
    text = AERO_DATA.read_text()
    assert re.search(r"ndust_spc\s*=\s*28\b", text)
    block = text[text.index("dust_spc( ndust_spc ) = (/") :]
    block = block[: block.index("/) ) /)") + len("/) ) /)")]
    rows = re.findall(
        r"emis_table\(\s*'[^']*',\s*'([^']*)',\s*([\d.Ee+-]+),\s*\(/\s*([\d.Ee+-]+),\s*([\d.Ee+-]+)\s*/\)",
        block,
    )
    source = tuple((name.strip(), float(mw), float(f), float(c)) for name, mw, f, c in rows)
    assert len(source) == 28
    copy = tuple((n, float(mw), float(f), float(c)) for n, mw, f, c in SOURCE_DUST_SPC)
    assert copy == source


def test_default_driver_splits_each_mode_into_the_28_source_species():
    """Without a caller table, silicon gets 19.425 % of fine dust and SOIL 95.358 % of coarse."""
    case = _port_case()
    case["settings"] = Ls99Settings()
    out = fengsha_dust_driver(**case)
    rates = out["species_rates"][0, 0]
    fine, coarse = out["mode_mass"][0, 0]
    assert rates.shape == (2, 28) and fine > 0 and coarse > 0
    si, soil = ls99.DUST_SPC_NAMES.index("SI"), ls99.DUST_SPC_NAMES.index("SOIL")
    assert rates[0, si] == fine * 0.19425
    assert rates[1, soil] == coarse * 0.95358
    assert rates[0, soil] == 0 and rates[1, si] == 0
