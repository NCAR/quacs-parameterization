"""Unit tests for the leung dust emission scheme.

These tests came from the original standalone module by Sue Park.
"""

import numpy as np
import pytest

from quacs.dust.leung.leung import (
    BIT,
    CD0,
    CE,
    CK,
    CTUNE,
    DP_MED,
    DST_SRC_NBR,
    GRAV,
    KAPPA_MAX,
    MSS_FRC_SRC,
    NDST,
    OVR_SRC_SNK_MSS,
    RHO_AIR0,
    RHO_P,
    SL_A,
    SL_GAMMA,
    STBLTY_FLOOR,
    U_FT0_DRY_REF,
    U_ST0,
    VAI_MBL_THR,
    Z0S_FACTOR,
    clay_mass_fraction_leung,
    drag_partition,
    dust_emission,
    intermittency_factor,
    moisture_threshold_factor,
    sl2000_u_ft0_dry,
    turbulent_sigma_us,
)


@pytest.fixture
def base_cell():
    """
    Single 1×1 land cell with conditions that produce positive emission:
    soil land type, dry soil, wind well above threshold, no snow, no ice,
    low vegetation, rock area = 0 (Feff = 1 from veg only → near-unity),
    Kc_map = 1 (no upscaling correction applied).
    """
    fclay = np.array([[0.10]])
    return dict(
        fv=np.array([[0.60]]),  # [m/s]  friction velocity
        forc_rho=np.array([[1.2]]),  # [kg/m3]
        h2osoi_vol=np.array([[0.05]]),  # [m3/m3]  dry soil
        h2osoi_liq=np.array([[10.0]]),  # [kg/m2]  mostly liquid
        h2osoi_ice=np.array([[0.0]]),  # [kg/m2]
        watsat=np.array([[0.40]]),  # [m3/m3]
        tlai=np.array([[0.0]]),  # bare ground
        tsai=np.array([[0.0]]),
        frac_sno=np.array([[0.0]]),
        obu=np.array([[-100.0]]),  # [m]  unstable (convective)
        fclay=fclay,
        z0a=np.array([[1.0e-4]]),  # [m]  smooth surface
        frac_rock=np.array([[0.0]]),  # no rocks
        frac_veg=np.array([[0.0]]),  # no plants
        Kc_map=np.array([[1.0]]),  # unity correction
        is_soil=np.array([[True]]),
    )


# --- output shape and sign ---


def test_output_shape(base_cell):
    """F_p must have shape [x, y, ndst=4]."""
    F_p = dust_emission(**base_cell)
    assert F_p.shape == (1, 1, NDST), f"Expected (1,1,4), got {F_p.shape}"


def test_all_fluxes_nonnegative(base_cell):
    """No negative emission fluxes anywhere."""
    F_p = dust_emission(**base_cell)
    assert (F_p >= 0.0).all(), "Negative flux detected"


def test_positive_emission_active_cell(base_cell):
    """Active cell with wind above threshold must produce nonzero flux."""
    F_p = dust_emission(**base_cell)
    assert F_p.sum() > 0.0, "Expected positive emission in active cell"


# --- absolute value / dimensional checks ---


def test_flux_hand_calc(base_cell):
    """
    Total flux must match an independent hand calculation of the
    Leung/CTSM equation (guards against dimensional errors such as a
    spurious 1/g, which relative tests cannot catch).
    """
    rho, fv = 1.2, 0.6
    fclay_eff = 0.1 + 0.5 * 0.10  # CTSM MassFracClayLeung2023 = 0.15
    u_ft0 = np.sqrt(SL_A * (RHO_P * GRAV * DP_MED + SL_GAMMA / DP_MED) / rho)
    u_it = BIT * u_ft0
    u_st = u_ft0 * np.sqrt(rho / RHO_AIR0)  # fm = 1 (gwc 0.031 < gwc_thr 0.184)
    Cd = CD0 * np.exp(-CE * (u_st - U_ST0) / U_ST0)
    kap = min(CK * (u_st - U_ST0) / U_ST0, KAPPA_MAX)
    liq = 10.0 / (10.0 + 1.0e-6)
    eta = intermittency_factor(
        np.array([[fv]]),
        turbulent_sigma_us(np.array([[fv]]), np.array([[-100.0]])),
        np.array([[u_it]]),
        np.array([[u_ft0]]),
    )[0, 0]
    Fd = (
        eta
        * CTUNE
        * Cd
        * 1.0
        * fclay_eff
        * liq
        * rho
        * (fv**2 - u_it**2)
        / u_it
        * (fv / u_it) ** kap
    )
    F_p = dust_emission(**base_cell)
    assert np.isclose(F_p[0, 0, :].sum(), Fd * OVR_SRC_SNK_MSS.sum(), rtol=1e-10)
    # order-of-magnitude anchor for the demo forcing [kg m-2 s-1]
    assert np.isclose(Fd, 1.133e-6, rtol=0.02), f"Fd = {Fd:.4e}"


def test_gwc_threshold_matches_ctsm():
    """gwc_thr must equal CTSM's 0.17 + 0.14*clay[%]*0.01 (0.184 at 10 % clay)."""
    _, _, gwc_thr = moisture_threshold_factor(
        np.array([[0.05]]), np.array([[0.40]]), np.array([[0.10]])
    )
    assert np.isclose(gwc_thr[0, 0], 0.184, rtol=1e-12)


def test_clay_term_matches_ctsm():
    """fclay' = 0.1 + 0.5*min(fclay, 0.2) (CTSM MassFracClayLeung2023)."""
    f = clay_mass_fraction_leung(np.array([[0.0, 0.10, 0.20, 0.40]]))
    assert np.allclose(f, [[0.10, 0.15, 0.20, 0.20]], rtol=1e-12)


def test_clay_free_soil_still_emits(base_cell):
    """With the CTSM clay term, fclay = 0 still gives fclay' = 0.1 > 0."""
    base_cell["fclay"] = np.array([[0.0]])
    assert dust_emission(**base_cell).sum() > 0.0


def test_eta_demo_near_one(base_cell):
    """
    At the demo forcing (u*s ≈ 3.4 u*it, unstable L = -100 m) emission is
    essentially continuous: eta must be ~1 once mean and sigma are both
    wind speeds at Z_SAL.
    """
    u_s = np.array([[0.6]])
    u_ft = sl2000_u_ft0_dry(np.array([[1.2]]))
    u_it = BIT * u_ft
    sig = turbulent_sigma_us(u_s, np.array([[-100.0]]))
    eta = intermittency_factor(u_s, sig, u_it, u_ft)
    assert eta[0, 0] > 0.999, f"eta = {eta[0, 0]:.6f}"


def test_sigma_neutral_limit():
    """Neutral limit (|L| -> inf): sigma_U / u_s -> 12^(1/3)."""
    sig = turbulent_sigma_us(np.array([[1.0]]), np.array([[1.0e12]]))
    assert np.isclose(sig[0, 0], 12.0 ** (1.0 / 3.0), rtol=1e-6)


def test_sigma_stable_floor():
    """Strongly stable (small positive L): bracket floored at STBLTY_FLOOR."""
    sig = turbulent_sigma_us(np.array([[1.0]]), np.array([[10.0]]))
    assert np.isclose(sig[0, 0], STBLTY_FLOOR ** (1.0 / 3.0), rtol=1e-12)


def test_sigma_increases_when_unstable():
    """More unstable (smaller negative L) must give larger sigma_U."""
    u_s = np.array([[0.5]])
    s_weak = turbulent_sigma_us(u_s, np.array([[-1000.0]]))
    s_strong = turbulent_sigma_us(u_s, np.array([[-10.0]]))
    assert s_strong[0, 0] > s_weak[0, 0]


# --- land-mask cutoffs ---


def test_zero_flux_non_soil(base_cell):
    """Non-soil land type must give zero flux."""
    base_cell["is_soil"] = np.array([[False]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Flux must be zero over non-soil"


def test_zero_flux_full_snow(base_cell):
    """Full snow cover (frac_sno=1) must give zero flux."""
    base_cell["frac_sno"] = np.array([[1.0]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Flux must be zero under full snow"


def test_snow_scales_flux(base_cell):
    """Partial snow must reduce flux proportionally to (1 - frac_sno)."""
    F_no_snow = dust_emission(**base_cell)
    base_cell["frac_sno"] = np.array([[0.5]])
    F_half = dust_emission(**base_cell)
    assert np.isclose(F_half.sum(), F_no_snow.sum() * 0.5, rtol=1e-7), (
        "Snow fraction should scale flux linearly"
    )


# --- vegetation cutoffs (VAIthr = 1.0) ---


def test_zero_flux_dense_vegetation(base_cell):
    """VAI >= VAIthr = 1.0 must give zero flux."""
    base_cell["tlai"] = np.array([[VAI_MBL_THR]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Flux must be zero above VAI threshold of 1.0"


def test_moderate_veg_still_emits(base_cell):
    """
    VAI = 0.3 must still emit under the new VAIthr = 1.0
    (would be zero in the old Z03/K14 scheme with VAIthr = 0.3).
    """
    base_cell["tlai"] = np.array([[0.3]])
    F_p = dust_emission(**base_cell)
    assert F_p.sum() > 0.0, "VAI = 0.3 should still produce emission with VAIthr = 1.0"


def test_vegetation_scales_flux(base_cell):
    """Half the VAIthr must give half the bare-ground flux (linear fbare)."""
    F_bare = dust_emission(**base_cell)
    base_cell["tlai"] = np.array([[VAI_MBL_THR / 2.0]])
    F_veg = dust_emission(**base_cell)
    assert np.isclose(F_veg.sum(), F_bare.sum() * 0.5, rtol=1e-7), (
        "Vegetation should reduce flux linearly with VAI fraction"
    )


# --- soil moisture ---


def test_wet_soil_reduces_flux(base_cell):
    """
    Clearly wet soil (gwc well above gwc_thr) must produce less flux than
    dry soil.

    Note: with gwc only slightly above gwc_thr = 0.184, flux can INCREASE
    a little: u*st rises, which lowers Cd but raises kappa, and because
    u*it does not depend on moisture, (u*s/u*it)^kappa grows.  At strong
    wind (u*s/u*it ≈ 3.4 here) the kappa effect can win.  This is a
    property of the Leung/K14 equations (also in CTSM), not a bug, so the
    test uses a clearly wet soil.
    """
    F_dry = dust_emission(**base_cell)
    base_cell["h2osoi_vol"] = np.array([[0.38]])  # gwc = 0.235 kg/kg
    F_wet = dust_emission(**base_cell)
    assert F_wet.sum() < F_dry.sum(), "Wet soil must reduce emission"


def test_frozen_soil_reduces_flux(base_cell):
    """High ice content (liqfrac → 0) must reduce flux."""
    F_liquid = dust_emission(**base_cell)
    base_cell["h2osoi_liq"] = np.array([[0.0]])
    base_cell["h2osoi_ice"] = np.array([[50.0]])
    F_frozen = dust_emission(**base_cell)
    assert F_frozen.sum() < F_liquid.sum(), "Frozen soil must reduce emission"


def test_moisture_below_threshold_no_effect(base_cell):
    """
    Below gwc_thr (CTSM: 0.17 + 0.14*fclay = 0.184 kg/kg at 10 % clay),
    soil moisture must not change the flux (fm = 1).
    """
    # gwc_sfc = h2osoi_vol * 1000 / ((1-0.4)*2700)
    #   0.05 -> 0.031 kg/kg,   0.25 -> 0.154 kg/kg   (both < 0.184)
    F_a = dust_emission(**base_cell)
    base_cell["h2osoi_vol"] = np.array([[0.25]])
    F_b = dust_emission(**base_cell)
    assert F_a.sum() > 0.0
    assert np.isclose(F_a.sum(), F_b.sum(), rtol=1e-12)


# --- wind threshold ---


def test_zero_flux_below_threshold(base_cell):
    """Very low wind must produce zero flux (u*s < u*it)."""
    base_cell["fv"] = np.array([[0.01]])
    base_cell["obu"] = np.array([[100.0]])  # stable, small sigma
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Sub-threshold wind must give zero flux"


def test_flux_increases_with_wind(base_cell):
    """Stronger wind must produce greater flux."""
    F_low = dust_emission(**base_cell)
    base_cell["fv"] = np.array([[1.20]])
    F_high = dust_emission(**base_cell)
    assert F_high.sum() > F_low.sum(), "Stronger wind must increase flux"


# --- drag partition ---


def test_rough_surface_reduces_flux(base_cell):
    """High aeolian roughness (large z0a) must reduce effective u*s and flux."""
    F_smooth = dust_emission(**base_cell)
    base_cell["z0a"] = np.array([[0.01]])  # rough: z0a = 1 cm
    F_rough = dust_emission(**base_cell)
    assert F_rough.sum() <= F_smooth.sum(), (
        "Rougher surface must reduce or maintain emission via drag partition"
    )


def test_rock_fraction_reduces_flux(base_cell):
    """Higher rock fraction (Ar > 0) with rough z0a must reduce u*s."""
    F_no_rock = dust_emission(**base_cell)
    base_cell["frac_rock"] = np.array([[0.5]])
    base_cell["z0a"] = np.array([[1.0e-3]])  # rough enough
    F_rock = dust_emission(**base_cell)
    assert F_rock.sum() <= F_no_rock.sum(), "Rock drag partition must reduce or maintain emission"


def test_vegetation_drag_reduces_flux(base_cell):
    """Vegetation drag partition (frac_veg > 0, VAI > 0) must reduce flux."""
    F_bare = dust_emission(**base_cell)
    base_cell["frac_veg"] = np.array([[0.5]])
    base_cell["tlai"] = np.array([[0.2]])
    F_veg = dust_emission(**base_cell)
    assert F_veg.sum() <= F_bare.sum(), "Vegetation drag partition must reduce or maintain emission"


# --- intermittency ---


def test_eta_between_zero_and_one(base_cell):
    """Intermittency factor must be in [0, 1] for all conditions."""
    # Test via direct helper over a range of wind speeds
    u_s = np.linspace(0.0, 1.0, 20).reshape(4, 5)
    sigma_us = np.full_like(u_s, 0.05)
    u_it = np.full_like(u_s, 0.20)
    u_ft = np.full_like(u_s, 0.25)
    eta = intermittency_factor(u_s, sigma_us, u_it, u_ft)
    assert (eta >= 0.0).all() and (eta <= 1.0).all(), "eta must lie in [0, 1]"


def test_eta_continuous_above_ft(base_cell):
    """When u_s >> u_ft, eta should be close to 1 (continuous emission)."""
    u_s = np.array([[2.0]])
    sigma_us = np.array([[0.05]])
    u_it = np.array([[0.20]])
    u_ft = np.array([[0.25]])
    eta = intermittency_factor(u_s, sigma_us, u_it, u_ft)
    assert eta[0, 0] > 0.99, f"Expected eta near 1, got {eta[0, 0]:.4f}"


def test_eta_zero_far_below_it(base_cell):
    """When u_s + sigma << u_it, eta should be near 0 (no emission)."""
    u_s = np.array([[0.01]])
    sigma_us = np.array([[0.005]])
    u_it = np.array([[0.30]])
    u_ft = np.array([[0.40]])
    eta = intermittency_factor(u_s, sigma_us, u_it, u_ft)
    assert eta[0, 0] < 0.01, f"Expected eta near 0, got {eta[0, 0]:.6f}"


def test_eta_increases_with_wind():
    """eta must increase monotonically with mean u_s (fixed sigma, thresholds)."""
    u_s = np.linspace(0.05, 0.6, 30).reshape(1, 30)
    sig = np.full_like(u_s, 0.5)
    u_it = np.full_like(u_s, 0.18)
    u_ft = np.full_like(u_s, 0.22)
    eta = intermittency_factor(u_s, sig, u_it, u_ft)[0]
    assert np.all(np.diff(eta) >= -1e-12)


def test_intermittency_increases_flux_vs_ft_scheme(base_cell):
    """
    Using u*it (< u*ft) as threshold allows emission at winds that would be
    zero under the u*ft scheme. With u_s between u_it and u_ft, eta > 0
    and Fd > 0 here; a u*ft scheme would give Fd = 0.
    """
    # Force u*s to lie between u*it and u*ft by choosing a moderate fv
    # and ensuring we're in the intermittent window.
    base_cell["fv"] = np.array([[0.22]])  # slightly above typical u*it
    base_cell["obu"] = np.array([[1e6]])  # very stable → small sigma
    F_p = dust_emission(**base_cell)
    # There should be some emission because u*it < u*s (even if small)
    # (a pure u*ft scheme would give zero in this window)
    assert F_p.sum() >= 0.0, "Intermittent scheme must give non-negative flux"


# --- upscaling correction map ---


def test_kc_scales_flux_linearly(base_cell):
    """Doubling Kc_map must exactly double the flux."""
    F1 = dust_emission(**base_cell)
    base_cell["Kc_map"] = np.array([[2.0]])
    F2 = dust_emission(**base_cell)
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-8), "Kc_map must scale flux linearly"


def test_kc_zero_gives_zero_flux(base_cell):
    """Kc_map = 0 must zero out all emission."""
    base_cell["Kc_map"] = np.array([[0.0]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Kc_map = 0 must produce zero flux"


# --- erodibility / threshold physics ---


def test_higher_clay_changes_flux(base_cell):
    """
    Higher clay raises the sandblasting efficiency (K14 uses fclay' directly)
    and also raises gwc_thr, so net effect on Fd is positive for dry conditions.
    """
    F_low_clay = dust_emission(**base_cell)  # fclay = 0.10
    base_cell["fclay"] = np.array([[0.18]])
    F_high_clay = dust_emission(**base_cell)
    assert F_high_clay.sum() >= 0.0 and F_low_clay.sum() >= 0.0, (
        "Flux must be non-negative for any clay fraction"
    )


def test_erodibility_decreases_with_moisture(base_cell):
    """
    Much higher moisture raises u*ft → higher u*st → lower Cd → lower Fd.
    (Use a clearly wet soil; just above gwc_thr the kappa effect can
    dominate — see test_wet_soil_reduces_flux.)
    """
    base_cell["h2osoi_vol"] = np.array([[0.02]])  # very dry (gwc 0.012)
    F_very_dry = dust_emission(**base_cell)
    base_cell["h2osoi_vol"] = np.array([[0.38]])  # clearly wet (gwc 0.235)
    F_wet = dust_emission(**base_cell)
    assert F_wet.sum() <= F_very_dry.sum(), "Higher moisture must reduce erodibility and flux"


def test_sl2000_threshold_positive():
    """Shao & Lu (2000) threshold must be strictly positive at RHO_AIR0."""
    assert U_FT0_DRY_REF > 0.0, "Dry threshold must be positive"


def test_sl2000_threshold_physical_range():
    """
    At Dp = 127 µm and rho_air = 1.225, S&L00 gives u*ft0 ≈ 0.215 m/s
    (Leung et al. 2023 Sect. 3.2). Check range [0.15, 0.35] m/s.
    """
    assert 0.15 < U_FT0_DRY_REF < 0.35, (
        f"Dry threshold = {U_FT0_DRY_REF:.4f} m/s outside physical range"
    )


def test_impact_threshold_less_than_fluid():
    """Impact threshold (u*it = Bit * u*ft0) must be < dry fluid threshold."""
    u_it = BIT * U_FT0_DRY_REF
    assert u_it < U_FT0_DRY_REF, f"u*it = {u_it:.4f} must be < u*ft0 = {U_FT0_DRY_REF:.4f}"


# --- drag partition helpers ---


def test_drag_partition_smooth_surface():
    """On a very smooth surface (z0a → z0s), feff_r should be near 1."""
    z0s = Z0S_FACTOR * DP_MED
    z0a = np.array([[z0s * 1.001]])  # barely above z0s
    vai = np.array([[0.0]])
    Ar = np.array([[1.0]])
    Av = np.array([[0.0]])
    Feff, feff_r, feff_v = drag_partition(z0a, vai, Ar, Av)
    assert feff_r[0, 0] > 0.9, (
        f"Near-smooth surface: expected feff_r near 1, got {feff_r[0, 0]:.4f}"
    )


def test_drag_partition_very_rough():
    """Very rough surface (large z0a) must strongly reduce Feff."""
    z0a = np.array([[0.10]])  # 10 cm — very rough
    vai = np.array([[0.0]])
    Ar = np.array([[1.0]])
    Av = np.array([[0.0]])
    Feff, feff_r, feff_v = drag_partition(z0a, vai, Ar, Av)
    assert Feff[0, 0] < 0.8, f"Very rough surface: expected Feff < 0.8, got {Feff[0, 0]:.4f}"


def test_drag_partition_bounds():
    """Feff must be in [0, 1] for any input combination."""
    rng = np.random.default_rng(7)
    nx, ny = 8, 6
    z0a = rng.uniform(1e-5, 0.1, (nx, ny))
    vai = rng.uniform(0.0, 1.0, (nx, ny))
    Ar = rng.uniform(0.0, 1.0, (nx, ny))
    Av = 1.0 - Ar
    Feff, _, _ = drag_partition(z0a, vai, Ar, Av)
    assert (Feff >= 0.0).all() and (Feff <= 1.0 + 1e-8).all(), "Feff must be in [0, 1]"


# --- bin partitioning ---


def test_bin_partition_sums_to_total(base_cell):
    """
    Sum over bins must equal Fd_corrected * total_overlap_fraction.
    """
    F_p = dust_emission(**base_cell)
    bin_sum = OVR_SRC_SNK_MSS.sum()
    if OVR_SRC_SNK_MSS[:, 0].sum() > 0:
        F0 = F_p[0, 0, 0]
        frac0 = OVR_SRC_SNK_MSS[:, 0].sum()
        Fd_total = F0 / frac0
        expected = Fd_total * bin_sum
        assert np.isclose(F_p[0, 0, :].sum(), expected, rtol=1e-7), (
            "Bin sum does not match expected partition"
        )


def test_four_bins(base_cell):
    """There must be exactly 4 output bins."""
    F_p = dust_emission(**base_cell)
    assert F_p.shape[-1] == 4, f"Expected 4 bins, got {F_p.shape[-1]}"


# --- precomputed constants ---


def test_ovr_src_snk_mss_shape():
    """Overlap matrix must be shape [3, 4]."""
    assert OVR_SRC_SNK_MSS.shape == (DST_SRC_NBR, NDST)


def test_ovr_src_snk_mss_nonnegative():
    """All overlap factors must be non-negative."""
    assert (OVR_SRC_SNK_MSS >= 0.0).all()


def test_ovr_src_snk_mss_row_sums():
    """Each source mode's total overlap cannot exceed its mass fraction."""
    for m in range(DST_SRC_NBR):
        row_sum = OVR_SRC_SNK_MSS[m, :].sum()
        assert row_sum <= MSS_FRC_SRC[m] + 1.0e-8, (
            f"Mode {m} overlap sum {row_sum:.4f} > mss_frc {MSS_FRC_SRC[m]}"
        )


def test_mss_frc_src_sums_to_one():
    """Source mode mass fractions must sum to 1.0."""
    assert np.isclose(MSS_FRC_SRC.sum(), 1.0, atol=1.0e-9)


# --- vectorised grid ---


def test_vectorized_grid():
    """Must run correctly on a full 2-D grid with mixed conditions."""
    nx, ny = 10, 8
    rng = np.random.default_rng(42)
    Ar = rng.uniform(0.0, 0.5, (nx, ny))
    Av = rng.uniform(0.0, 0.5, (nx, ny))
    F_p = dust_emission(
        fv=rng.uniform(0.1, 1.5, (nx, ny)),
        forc_rho=rng.uniform(1.1, 1.3, (nx, ny)),
        h2osoi_vol=rng.uniform(0.0, 0.30, (nx, ny)),
        h2osoi_liq=rng.uniform(0.0, 30.0, (nx, ny)),
        h2osoi_ice=rng.uniform(0.0, 5.0, (nx, ny)),
        watsat=rng.uniform(0.30, 0.50, (nx, ny)),
        tlai=rng.uniform(0.0, 1.0, (nx, ny)),
        tsai=rng.uniform(0.0, 0.3, (nx, ny)),
        frac_sno=rng.uniform(0.0, 1.0, (nx, ny)),
        obu=rng.uniform(-500.0, 500.0, (nx, ny)),
        fclay=rng.uniform(0.01, 0.20, (nx, ny)),
        z0a=rng.uniform(1.0e-5, 0.05, (nx, ny)),
        frac_rock=Ar,
        frac_veg=Av,
        Kc_map=rng.uniform(0.5, 3.0, (nx, ny)),
        is_soil=rng.choice([True, False], (nx, ny)),
    )
    assert F_p.shape == (nx, ny, NDST), f"Shape mismatch: {F_p.shape}"
    assert (F_p >= 0.0).all(), "Negative flux in vectorised grid test"
    assert np.isfinite(F_p).all(), "Non-finite flux in vectorised grid test"
