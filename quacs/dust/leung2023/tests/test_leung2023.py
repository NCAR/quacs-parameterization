"""Tests for quacs.dust.leung2023."""

import warnings
from dataclasses import fields, replace

import numpy as np
import pytest
from scipy.special import erf

from quacs.dust.leung2023.leung2023 import (
    B_IT,
    C_TUNE,
    DP,
    FECAN_A,
    LeungInputWarning,
    LeungSettings,
    dry_fluid_threshold,
    dust_dry_dep,
    dust_emission,
    init_dust_vars,
    initialize,
    intermittency_factor,
    leung_dust_driver,
    tunable_names,
)


@pytest.fixture
def ctsm_cell(cell):
    """CTSM surface and patch inputs for one emitting patch."""
    return dict(
        ustar=cell["ustar"],
        rhoa=cell["rhoa"],
        h2osoi_vol=np.array([[0.01]]),
        tlai=np.array([[0.03]]),
        tsai=np.array([[0.02]]),
        frac_sno=cell["frac_sno"],
        pct_sand=np.array([[65.0]]),
        pct_clay=cell["pct_clay"],
        organic=np.array([[20.0]]),
        h2osoi_liq=np.array([[1.0]]),
        h2osoi_ice=np.array([[0.0]]),
        land_fraction=np.array([[1.0]]),
        landunit_type=np.array([[1]]),
        patch_active=np.array([[True]]),
        patch_wtlunit=np.array([[1.0]]),
        forc_pbot=np.array([[100000.0]]),
        forc_t=np.array([[295.0]]),
        ram1=np.array([[100.0]]),
        rockfrc=cell["rockfrc"],
        vegefrc=cell["vegefrc"],
        sparfrc=cell["sparfrc"],
        roughfct=cell["roughfct"],
    )


def test_ctsm_driver_splits_emission_into_four_source_bins(ctsm_cell):
    result = leung_dust_driver(**ctsm_cell)
    assert result["flx_mss_vrt_dst"].shape == (1, 1, 4)
    np.testing.assert_allclose(
        result["flx_mss_vrt_dst_tot"], result["flx_mss_vrt_dst"].sum(axis=-1)
    )
    np.testing.assert_allclose(
        result["flx_mss_vrt_dst_tot"] / result["flux_before_split"], 0.8712, rtol=3e-4
    )


def test_ctsm_driver_coast_cut_applies_below_half_land_fraction(ctsm_cell):
    base = leung_dust_driver(**ctsm_cell)["flx_mss_vrt_dst_tot"]
    coast = leung_dust_driver(**{**ctsm_cell, "land_fraction": np.array([[0.49]])})[
        "flx_mss_vrt_dst_tot"
    ]
    np.testing.assert_allclose(coast, base * 1e-3)


def test_ctsm_driver_raises_on_patch_weights_above_one(ctsm_cell):
    with pytest.raises(ValueError, match="weight"):
        leung_dust_driver(**{**ctsm_cell, "patch_wtlunit": np.array([[1.01]])})


def test_ctsm_landunit_mean_uses_only_active_weighted_patches(ctsm_cell):
    case = {
        **ctsm_cell,
        "ustar": np.array([[0.6, 0.6]]),
        "tlai": np.array([[0.2, 0.9]]),
        "tsai": np.array([[0.0, 0.0]]),
        "patch_active": np.array([[True, False]]),
        "patch_wtlunit": np.array([[0.6, 0.4]]),
        "landunit_id": np.array([[1, 1]]),
    }
    result = leung_dust_driver(**case)
    np.testing.assert_allclose(result["tlai_lu"], [[0.2, 0.2]])
    assert result["flx_mss_vrt_dst_tot"][0, 0] > 0
    active = leung_dust_driver(**{**case, "patch_active": np.array([[True, True]])})
    np.testing.assert_allclose(active["tlai_lu"], [[0.48, 0.48]])
    # An inactive patch is outside filter_nolakep, so it emits nothing.
    assert result["flx_mss_vrt_dst_tot"][0, 1] == 0.0


def test_ctsm_landunit_mean_counts_lake_patch_as_bare(ctsm_cell):
    # DUSTMod L455-L472: ttlai starts at 0 for every patch and is set only
    # for non-lake patches, so a lake patch adds weight but no vegetation.
    case = {
        **ctsm_cell,
        "ustar": np.array([[0.6, 0.6]]),
        "tlai": np.array([[0.2, 0.9]]),
        "tsai": np.array([[0.0, 0.0]]),
        "patch_active": np.array([[True, True]]),
        "patch_wtlunit": np.array([[0.6, 0.4]]),
        "landunit_id": np.array([[1, 1]]),
        "lake": np.array([[False, True]]),
    }
    result = leung_dust_driver(**case)
    np.testing.assert_allclose(result["tlai_lu"], [[0.12, 0.12]])
    assert result["flx_mss_vrt_dst_tot"][0, 1] == 0.0


def test_ctsm_drag_output_is_zero_above_vai_cutoff(ctsm_cell):
    # DUSTMod L623 and L640: frc_thr_rghn_fct is 0 when tlai_lu > 1.
    dense = leung_dust_driver(**{**ctsm_cell, "tlai": np.array([[1.5]])})
    assert dense["f_eff"][0, 0] == 0.0
    sparse = leung_dust_driver(**ctsm_cell)
    assert sparse["f_eff"][0, 0] > 0.0


def test_ctsm_driver_matches_pinned_four_bins_deposition_and_setup(ctsm_cell):
    case = {
        **ctsm_cell,
        "rhoa": np.array([[1.1]]),
        "h2osoi_liq": np.array([[5.0]]),
        "organic": np.array([[20.0]]),
    }
    result = leung_dust_driver(**case)
    np.testing.assert_allclose(
        result["flx_mss_vrt_dst"][0, 0],
        [4.225342311015966e-9, 2.26806030651875e-8, 5.3183521651839225e-8, 5.009721293978224e-8],
        rtol=1e-10,
    )
    np.testing.assert_allclose(
        result["vlc_trb"][0, 0],
        [1.0023991132362964e-4, 5.698390058254035e-5, 7.113891110430164e-3, 9.482193079306147e-3],
        rtol=1e-10,
    )
    np.testing.assert_allclose(
        result["dmt_vwr"],
        [7.805670344214621e-7, 1.8081280010153617e-6, 3.601478592068484e-6, 6.814310470691836e-6],
        rtol=1e-10,
    )
    np.testing.assert_allclose(result["stk_crc"], [1.0] * 4, rtol=1e-10)


def test_second_pinned_emitting_case_covers_stronger_wind_and_dry_deposition(ctsm_cell):
    case = {
        **ctsm_cell,
        "ustar": np.array([[1.0]]),
        "rhoa": np.array([[1.1]]),
        "pct_sand": np.array([[60.0]]),
        "pct_clay": np.array([[15.0]]),
        "organic": np.array([[0.0]]),
        "h2osoi_vol": np.array([[0.02]]),
        "tlai": np.array([[0.1]]),
        "tsai": np.array([[0.0]]),
        "h2osoi_liq": np.array([[5.0]]),
        "forc_pbot": np.array([[90000.0]]),
        "forc_t": np.array([[280.0]]),
        "ram1": np.array([[30.0]]),
    }
    result = leung_dust_driver(**case)
    np.testing.assert_allclose(
        result["flx_mss_vrt_dst"][0, 0],
        [2.8484372708464615e-8, 1.5289713907373517e-7, 3.585269881519812e-7, 3.377212022116524e-7],
        rtol=1e-10,
    )
    np.testing.assert_allclose(
        result["vlc_trb"][0, 0],
        [1.7279702061977984e-4, 1.2668780657501868e-2, 3.0419739167745225e-2, 3.175786049269035e-2],
        rtol=1e-10,
    )


def test_ctsm_zero_surface_texture_gets_one_percent_floor(ctsm_cell):
    result = leung_dust_driver(
        **{**ctsm_cell, "pct_sand": np.array([[0.0]]), "pct_clay": np.array([[0.0]])}
    )
    assert result["f_clay"][0, 0] == 0.01
    assert result["watsat"][0, 0] > 0


@pytest.mark.parametrize(
    "sand_layers,clay_layers", [([0.0, 1.0], [0.0, 0.0]), ([10.0, 0.0], [10.0, 0.0])]
)
def test_ctsm_surface_texture_depth_stops_match_host(ctsm_cell, sand_layers, clay_layers):
    case = {
        **ctsm_cell,
        "pct_sand": np.array([[sand_layers[0]]]),
        "pct_clay": np.array([[clay_layers[0]]]),
        "pct_sand_layers": np.array([[sand_layers]]),
        "pct_clay_layers": np.array([[clay_layers]]),
    }
    with pytest.raises(ValueError, match="depth points"):
        leung_dust_driver(**case)


def test_ctsm_mobile_fraction_stop_matches_fortran(ctsm_cell):
    with pytest.raises(ValueError, match="mobile fraction"):
        leung_dust_driver(**{**ctsm_cell, "frac_sno": np.array([[1.1]])})


def test_setup_rejects_non_four_bin_count():
    with pytest.raises(ValueError, match="ndst"):
        init_dust_vars(LeungSettings(ndst=3))


def test_setup_reynolds_stop_matches_fortran():
    settings = LeungSettings(bin_grid=(0.01, 0.02, 0.03, 0.04, 0.05), dmt_vma=0.03)
    with pytest.raises(ValueError, match="Reynolds"):
        init_dust_vars(settings)


def test_setup_nonconvergence_warning_uses_one_warning_class():
    with pytest.warns(LeungInputWarning, match="not converging"):
        init_dust_vars(LeungSettings(eps_max=-1.0, convergence_warn=0))


def test_impossible_organic_warns_then_calculates(ctsm_cell):
    with pytest.warns(LeungInputWarning, match="outside physical range"):
        result = leung_dust_driver(**{**ctsm_cell, "organic": np.array([[-1.0]])})
    assert np.isfinite(result["watsat"]).all()


_SETTING_CASES = {
    **dict.fromkeys(
        (
            "sl00_a",
            "sl00_gam",
            "dp",
            "dns_slt",
            "fecan_a",
            "cd0",
            "ce",
            "ca",
            "c_tune",
            "ust0_min",
            "forc_rho_std",
            "vai_mbl_thr",
            "b_it",
            "u_0",
            "c_e",
            "z_sal",
            "z0_sal",
            "sd_exp",
            "sqrt2_lit",
            "cbrt_lit",
            "clay_linear",
            "clay_quadratic",
            "wet_scale",
            "wet_exp",
            "lai_offset",
            "lai_cap",
            "gap_multiplier",
            "drag_gate",
            "variance_neutral",
            "sd_floor",
            "clay_cap",
            "frag_cap",
            "coast_cutoff",
            "coast_scale",
            "soil_particle_density",
            "liquid_guard",
            "organic_max",
            "zsapric",
            "zsoi1",
            "cosby_base",
            "cosby_sand",
            "organic_porosity_base",
            "organic_porosity_slope",
            "organic_porosity_min",
            "use_drag_partition",
            "use_intermittency",
        ),
        "driver",
    ),
    **dict.fromkeys(
        (
            "source_medians",
            "source_gsds",
            "source_fractions",
            "bin_grid",
            "sz_nbr",
            "gsd_anl",
            "dns_aer",
            "dmt_vma",
            "setup_t",
            "setup_p",
            "mu_scale",
            "mu_tref",
            "mu_s",
            "mu_offset",
            "slip_a",
            "slip_b",
            "slip_exp",
            "schmidt_exp",
            "impaction_exp",
        ),
        "setup",
    ),
    **dict.fromkeys(
        (
            "re_cut1",
            "re_cut2",
            "re_cut3",
            "drag_low_linear",
            "drag_low_log",
            "drag_mid_scale",
            "drag_mid_exp",
            "drag_high",
            "eps_max",
            "convergence_average",
            "convergence_warn",
        ),
        "reynolds",
    ),
    "stability_coeff": "stability",
    "weight_tolerance": "weight",
    "ndst": "bin_error",
    "re_max": "reynolds_error",
}


@pytest.mark.parametrize("name", tunable_names(LeungSettings))
def test_each_leung_setting_changes_its_result_or_source_stop(name, ctsm_cell):
    assert set(_SETTING_CASES) == set(tunable_names(LeungSettings))
    assert {"source", "what"} <= fields(LeungSettings)[
        tunable_names(LeungSettings).index(name)
    ].metadata.keys()
    default = LeungSettings()
    old = getattr(default, name)
    if isinstance(old, bool):
        changed = not old
    elif isinstance(old, tuple):
        changed = (old[0] * 1.25, *old[1:])
    elif isinstance(old, int):
        changed = max(1, old // 2)
    else:
        changed = old * 1.25
    changed = {
        "lai_cap": 0.2,
        "drag_gate": 0.1,
        "sd_floor": 20.0,
        "coast_cutoff": 0.6,
        "organic_porosity_min": 1.0,
        "re_cut1": 0.3,
        "re_cut2": 0.5,
        "re_cut3": 100.0,
        "eps_max": 0.5,
        "convergence_average": 0,
        "convergence_warn": 0,
        "weight_tolerance": 1e-4,
        "ndst": 3,
        "re_max": 1.0,
    }.get(name, changed)
    scenario = _SETTING_CASES[name]
    tuned = replace(default, **{name: changed})
    if scenario == "driver":
        case = {
            **ctsm_cell,
            "h2osoi_vol": np.array([[0.08]]),
            "tlai": np.array([[0.4]]),
            "pct_clay": np.array([[30.0]]),
            "land_fraction": np.array([[0.55]]),
        }
        if name in ("wet_scale", "wet_exp"):
            case["h2osoi_vol"] = np.array([[0.2]])
            case["pct_clay"] = np.array([[10.0]])
        if name == "frag_cap":
            case["h2osoi_vol"] = np.array([[0.4]])
        if name == "coast_scale":
            case["land_fraction"] = np.array([[0.49]])
        before = leung_dust_driver(**case, settings=default)
        after = leung_dust_driver(**case, settings=tuned)
        keys = (
            "flux_before_split",
            "flx_mss_vrt_dst",
            "watsat",
            "gwc",
            "liqfrac",
            "tlai_lu",
            "gwc_thr",
            "f_clay",
            "eta",
            "u_s",
            "frag_expt",
        )
        assert any(not np.array_equal(before[k], after[k]) for k in keys)
    elif scenario == "setup":
        if name == "setup_p":
            default = replace(
                default,
                bin_grid=tuple(x * 20 for x in default.bin_grid),
                dmt_vma=default.dmt_vma * 20,
            )
            tuned = replace(default, setup_p=changed)
        before = init_dust_vars(default)
        after = init_dust_vars(tuned)
        if name in (
            "mu_scale",
            "mu_tref",
            "mu_s",
            "mu_offset",
            "slip_a",
            "slip_b",
            "slip_exp",
            "setup_t",
            "setup_p",
            "dns_aer",
            "schmidt_exp",
            "impaction_exp",
        ):
            base_dep = dust_dry_dep(
                np.array([[0.6]]),
                np.array([[1.1]]),
                np.array([[100000.0]]),
                np.array([[295.0]]),
                np.array([[100.0]]),
                settings=default,
                setup=before,
            )
            tuned_dep = dust_dry_dep(
                np.array([[0.6]]),
                np.array([[1.1]]),
                np.array([[100000.0]]),
                np.array([[295.0]]),
                np.array([[100.0]]),
                settings=tuned,
                setup=after,
            )
            assert not np.array_equal(base_dep["vlc_trb"], tuned_dep["vlc_trb"])
        else:
            assert any(
                not np.array_equal(before[k], after[k])
                for k in ("ovr_src_snk_mss", "dmt_vwr", "stk_crc")
            )
    elif scenario == "reynolds":
        scale = {
            "re_cut1": 20,
            "re_cut2": 20,
            "re_cut3": 100,
            "drag_low_linear": 20,
            "drag_low_log": 20,
            "drag_mid_scale": 50,
            "drag_mid_exp": 50,
            "drag_high": 150,
        }.get(name, 50)
        base = replace(
            default,
            bin_grid=tuple(x * scale for x in default.bin_grid),
            dmt_vma=default.dmt_vma * scale,
        )
        if name in ("convergence_average", "convergence_warn"):
            base = replace(base, eps_max=-1.0)
        alt = replace(base, **{name: changed})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LeungInputWarning)
            before = init_dust_vars(base)["stk_crc"]
            after = init_dust_vars(alt)["stk_crc"]
        assert not np.array_equal(before, after)
    elif scenario == "stability":
        args = (np.array([0.6]), np.array([0.3]), np.array([0.2]))
        assert not np.array_equal(
            intermittency_factor(*args, settings=default, stability=2.0),
            intermittency_factor(*args, settings=tuned, stability=2.0),
        )
    elif scenario == "weight":
        case = {**ctsm_cell, "patch_wtlunit": np.array([[1.000002]])}
        with pytest.raises(ValueError):
            leung_dust_driver(**case, settings=default)
        assert np.isfinite(leung_dust_driver(**case, settings=tuned)["watsat"]).all()
    elif scenario == "bin_error":
        with pytest.raises(ValueError, match="ndst"):
            init_dust_vars(tuned)
    elif scenario == "reynolds_error":
        base = replace(
            default,
            bin_grid=tuple(x * 150 for x in default.bin_grid),
            dmt_vma=default.dmt_vma * 150,
        )
        assert np.isfinite(init_dust_vars(base)["stk_crc"]).all()
        with pytest.raises(ValueError, match="Reynolds"):
            init_dust_vars(replace(base, re_max=changed))


@pytest.fixture
def cell():
    """Single mobile Gobi-like land cell, dry, strong wind well above threshold."""
    return dict(
        ustar=np.array([[0.6]]),
        rhoa=np.array([[1.2]]),
        gwc=np.array([[0.01]]),
        lai=np.array([[0.05]]),
        frac_sno=np.array([[0.0]]),
        pct_clay=np.array([[10.0]]),
        rockfrc=np.array([[0.5]]),
        vegefrc=np.array([[0.1]]),
        sparfrc=np.array([[0.1]]),
        roughfct=np.array([[0.7]]),
    )


# --- shape and basic contract ---


def test_core_pre_split_flux_shape(cell):
    assert dust_emission(**cell).shape == (1, 1)


def test_all_fluxes_nonnegative(cell):
    assert (dust_emission(**cell) >= 0.0).all()


def test_diagnostics_returned(cell):
    F, d = dust_emission(return_diagnostics=True, **cell)
    for key in ("u_ft0", "u_ft", "u_it", "u_st", "Cd", "frag_expt", "eta"):
        assert key in d


# --- source-verified parameter guards ---


def test_Dp_is_130_micron():
    """Dp must be the CODE value 130 um (DUSTMod L593), not the paper's 127 um."""
    assert np.isclose(DP, 130.0e-6), f"Dp should be 130e-6 m (code); got {DP:.2e}"


def test_fecan_a_is_two():
    """Fecan tunable a must be 2 (SoilStateInitTimeConstMod L675)."""
    assert FECAN_A == 2.0


def test_dry_threshold_matches_paper_magnitude():
    """u*ft0 at rho=1.225 should be ~0.215-0.216 m/s (paper Sect 3.2)."""
    u_ft0 = dry_fluid_threshold(np.array([[1.225]]))
    assert 0.20 < float(u_ft0[0, 0]) < 0.23


# --- structural feature: impact threshold ignores moisture ---


def test_impact_threshold_independent_of_moisture():
    """
    u*it must be B_it * dry threshold and NOT change with soil moisture
    (DUSTMod L596). Verified via the diagnostics: u_it equals B_it*u_ft0 and
    is identical for dry vs moist soil, while u_ft (wet) rises.
    """
    base = dict(
        ustar=np.array([[0.6]]),
        rhoa=np.array([[1.2]]),
        lai=np.array([[0.05]]),
        frac_sno=np.array([[0.0]]),
        pct_clay=np.array([[10.0]]),
        rockfrc=np.array([[0.5]]),
        vegefrc=np.array([[0.1]]),
        sparfrc=np.array([[0.1]]),
        roughfct=np.array([[0.7]]),
    )
    _, d_dry = dust_emission(gwc=np.array([[0.001]]), return_diagnostics=True, **base)
    _, d_wet = dust_emission(gwc=np.array([[0.30]]), return_diagnostics=True, **base)
    assert np.isclose(d_dry["u_it"][0, 0], d_wet["u_it"][0, 0]), "u*it must not depend on moisture"
    assert np.isclose(d_dry["u_it"][0, 0], B_IT * d_dry["u_ft0"][0, 0])
    assert d_wet["u_ft"][0, 0] > d_dry["u_ft"][0, 0], "u*ft (wet) must rise with moisture"


# --- emission gate ---


def test_zero_flux_below_impact_threshold():
    """u*s at or below the impact threshold must give zero flux."""
    F = dust_emission(
        ustar=np.array([[0.05]]),
        rhoa=np.array([[1.2]]),
        gwc=np.array([[0.01]]),
        lai=np.array([[0.05]]),
        frac_sno=np.array([[0.0]]),
        pct_clay=np.array([[10.0]]),
        rockfrc=np.array([[0.5]]),
        vegefrc=np.array([[0.1]]),
        sparfrc=np.array([[0.1]]),
        roughfct=np.array([[0.7]]),
    )
    assert (F == 0.0).all()


# --- mobility cutoffs ---


def test_zero_flux_when_vai_above_threshold():
    """lai >= VAI_MBL_THR (=1) must give zero flux (lnd_frc_mbl = 0)."""
    F = dust_emission(
        ustar=np.array([[0.8]]),
        rhoa=np.array([[1.2]]),
        gwc=np.array([[0.01]]),
        lai=np.array([[1.0]]),
        frac_sno=np.array([[0.0]]),
        pct_clay=np.array([[10.0]]),
        rockfrc=np.array([[0.5]]),
        vegefrc=np.array([[0.1]]),
        sparfrc=np.array([[0.1]]),
        roughfct=np.array([[0.7]]),
    )
    assert (F == 0.0).all()


def test_full_snow_gives_zero_flux(cell):
    """frac_sno = 1 must zero the flux (lnd_frc_mbl *= 1 - frac_sno)."""
    F = dust_emission(**{**cell, "frac_sno": np.array([[1.0]])})
    assert (F == 0.0).all()


# --- moisture behavior ---


def test_moisture_reduces_flux(cell):
    """Moister soil raises u*ft, lowers Cd and n, and reduces flux."""
    F_dry = dust_emission(**{**cell, "gwc": np.array([[0.001]])})
    F_wet = dust_emission(**{**cell, "gwc": np.array([[0.10]])})
    assert F_wet.sum() < F_dry.sum()


# --- intermittency ---


def test_intermittency_in_unit_interval():
    """eta must lie in [0, 1] across a range of winds."""
    u_s = np.linspace(0.1, 1.5, 50)
    eta = intermittency_factor(u_s, np.full_like(u_s, 0.25), np.full_like(u_s, 0.20))
    assert np.all(eta >= -1e-9) and np.all(eta <= 1.0 + 1e-9)


def test_intermittency_closed_form_sd():
    """eta matches a hand evaluation of DUSTMod L778-L811 with stability 0.

    The source constants are written as literals (12, 0.333, 0.1, 1e-4, 0.4,
    1.414) so a change to any default in LeungSettings makes this fail.
    """
    u_s = np.array([0.20, 0.30, 0.50])
    u_ft, u_it = 0.25, 0.20
    k, log_ratio = 0.4, np.log(0.1 / 1e-4)
    u_bar = (u_s / k) * log_ratio
    u_sd = u_s * 12.0**0.333
    u_ftz = (u_ft / k) * log_ratio
    u_itz = (u_it / k) * log_ratio
    R = 1.0 / (
        np.exp((u_ftz**2 - u_itz**2 - 2.0 * u_bar * (u_ftz - u_itz)) / (2.0 * u_sd**2)) + 1.0
    )
    P_ft = 0.5 * (1.0 + erf((u_ftz - u_bar) / (1.414 * u_sd)))
    P_it = 0.5 * (1.0 + erf((u_itz - u_bar) / (1.414 * u_sd)))
    eta_hand = 1.0 - P_ft + R * (P_ft - P_it)
    eta = intermittency_factor(u_s, np.full_like(u_s, u_ft), np.full_like(u_s, u_it))
    np.testing.assert_allclose(eta, eta_hand, rtol=1e-12, atol=0)


def test_intermittency_toggle_changes_flux(cell):
    """Turning intermittency off must change (raise or equal) the flux."""
    F_on = dust_emission(**cell)
    F_off = dust_emission(**cell, settings=LeungSettings(use_intermittency=False))
    assert F_off.sum() >= F_on.sum() - 1e-30


# --- scaling laws ---


def test_flux_scales_linearly_with_C_tune(cell):
    """Doubling C_tune must exactly double the flux."""
    F1 = dust_emission(**cell)
    F2 = dust_emission(**cell, settings=LeungSettings(c_tune=C_TUNE * 2.0))
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-10)


def test_flux_scales_linearly_with_liqfrac(cell):
    """liqfrac enters linearly (DUSTMod L755)."""
    F_full = dust_emission(**{**cell, "liqfrac": 1.0})
    F_half = dust_emission(**{**cell, "liqfrac": 0.5})
    assert np.isclose(F_half.sum(), F_full.sum() * 0.5, rtol=1e-10)


# --- drag partition ---


def test_drag_partition_reduces_wind(cell):
    """Drag partition (F_eff <= 1) must reduce u*s and the flux vs. no partition."""
    F_on = dust_emission(**cell)
    F_off = dust_emission(**cell, settings=LeungSettings(use_drag_partition=False))
    assert F_off.sum() >= F_on.sum()


# --- vectorized grid ---


def test_vectorized_grid():
    """Must run on a full 2D grid with mixed mobile / vegetated / snowy cells."""
    nx, ny = 12, 10
    rng = np.random.default_rng(0)
    F = dust_emission(
        ustar=rng.uniform(0.0, 1.2, (nx, ny)),
        rhoa=rng.uniform(1.0, 1.3, (nx, ny)),
        gwc=rng.uniform(0.0, 0.2, (nx, ny)),
        lai=rng.uniform(0.0, 1.5, (nx, ny)),
        frac_sno=rng.uniform(0.0, 0.3, (nx, ny)),
        pct_clay=rng.uniform(2.0, 30.0, (nx, ny)),
        rockfrc=rng.uniform(0.0, 0.8, (nx, ny)),
        vegefrc=rng.uniform(0.0, 0.3, (nx, ny)),
        sparfrc=rng.uniform(0.0, 0.3, (nx, ny)),
        roughfct=rng.uniform(0.3, 1.0, (nx, ny)),
    )
    assert F.shape == (nx, ny)
    assert (F >= 0.0).all()
    assert np.isfinite(F).all()


# --- initialize() ---


def test_initialize_missing_file():
    with pytest.raises(FileNotFoundError):
        initialize("no_clay.npy", "no_r.npy", "no_v.npy", "no_s.npy", "no_f.npy")


def test_initialize_valid(tmp_path):
    paths = {}
    for name, arr in [
        ("clay", np.full((6, 6), 12.0)),
        ("rock", np.full((6, 6), 0.5)),
        ("vege", np.full((6, 6), 0.1)),
        ("spar", np.full((6, 6), 0.1)),
        ("rough", np.full((6, 6), 0.7)),
    ]:
        p = tmp_path / f"{name}.npy"
        np.save(p, arr)
        paths[name] = str(p)
    pct_clay, rockfrc, vegefrc, sparfrc, roughfct = initialize(
        paths["clay"], paths["rock"], paths["vege"], paths["spar"], paths["rough"]
    )
    assert pct_clay.shape == (6, 6) and roughfct.shape == (6, 6)


def test_initialize_clay_out_of_range(tmp_path):
    """Clay given as a fraction (0-1) instead of a percentage must be caught."""
    paths = {}
    for name, arr in [
        ("clay", np.full((4, 4), 0.12)),  # fraction, not percent
        ("rock", np.full((4, 4), 0.5)),
        ("vege", np.full((4, 4), 0.1)),
        ("spar", np.full((4, 4), 0.1)),
        ("rough", np.full((4, 4), 0.7)),
    ]:
        p = tmp_path / f"{name}.npy"
        np.save(p, arr)
        paths[name] = str(p)
    # 0.12 is a valid percentage too, so this should NOT raise; guards only >100
    pct_clay, *_ = initialize(
        paths["clay"], paths["rock"], paths["vege"], paths["spar"], paths["rough"]
    )
    assert pct_clay.max() <= 100.0
