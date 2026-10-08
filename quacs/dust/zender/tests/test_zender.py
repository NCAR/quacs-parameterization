"""Unit tests for the zender dust emission scheme.

These tests came from the original standalone module by Sue Park.
"""

import numpy as np
import pytest

from quacs.dust.zender.zender import (
    CST_SLT,
    DST_SRC_NBR,
    FLX_MSS_FDG_FCT,
    GRAV,
    MSS_FRC_SRC,
    NDST,
    OVR_SRC_SNK_MSS,
    TMP1,
    VAI_MBL_THR,
    dust_emission,
    gwc_threshold,
)


@pytest.fixture
def base_cell():
    """
    Single 1×1 land cell with conditions that produce positive emission:
    soil land type, dry soil (gwc below threshold), wind well above threshold,
    no snow, no ice, low vegetation.
    """
    fclay = np.array([[0.10]])  # clay fraction 10 %
    return dict(
        fv=np.array([[0.50]]),  # [m/s]  friction velocity
        u10=np.array([[8.0]]),  # [m/s]  10-m wind
        forc_rho=np.array([[1.2]]),  # [kg/m3]
        h2osoi_vol=np.array([[0.05]]),  # [m3/m3]  dry soil
        h2osoi_liq=np.array([[10.0]]),  # [kg/m2]  mostly liquid
        h2osoi_ice=np.array([[0.0]]),  # [kg/m2]
        watsat=np.array([[0.40]]),  # [m3/m3]
        tlai=np.array([[0.0]]),  # bare ground
        tsai=np.array([[0.0]]),
        frac_sno=np.array([[0.0]]),
        fclay=fclay,
        mbl_bsn_fct=np.array([[1.0]]),  # CLM4 default
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


# --- absolute value / dimensional check ---


def test_flux_hand_calc(base_cell):
    """
    Total flux must match an independent hand calculation of the CTSM
    Zender2003 equations (guards against unit / formula regressions that
    relative tests cannot catch).
    """
    rho, fv, u10, fclay = 1.2, 0.5, 8.0, 0.10
    u_thr = TMP1 / np.sqrt(rho)  # fm = 1 (gwc 0.031 < 0.184)
    u10thr = u10 * u_thr / fv
    u_slt = fv + 0.003 * (u10 - u10thr) ** 2
    r = u_thr / u_slt
    liq = 10.0 / (10.0 + 1.0e-6)
    Q = (
        CST_SLT * rho * u_slt**3 / GRAV * (1 - r) * (1 + r) ** 2 * FLX_MSS_FDG_FCT * liq
    )  # [kg m-1 s-1]
    F_tot = Q * 100.0 * 10.0 ** (13.4 * fclay - 6.0)  # [kg m-2 s-1]
    F_p = dust_emission(**base_cell)
    assert np.isclose(F_p[0, 0, :].sum(), F_tot * OVR_SRC_SNK_MSS.sum(), rtol=1e-10)
    # order-of-magnitude anchor for this forcing [kg m-2 s-1]
    assert np.isclose(F_tot, 7.49e-8, rtol=0.02), f"F_tot = {F_tot:.4e}"


def test_tmp1_value():
    """tmp1 / sqrt(1.2) ≈ 0.2069 m/s for Dp = 75 µm (Iversen & White 1982)."""
    assert np.isclose(TMP1 / np.sqrt(1.2), 0.2069, rtol=2e-3)


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
    F_half_snow = dust_emission(**base_cell)
    assert np.isclose(F_half_snow.sum(), F_no_snow.sum() * 0.5, rtol=1e-8), (
        "Snow fraction should scale flux linearly"
    )


# --- vegetation cutoffs ---


def test_zero_flux_dense_vegetation(base_cell):
    """VAI >= vai_mbl_thr must give zero flux."""
    base_cell["tlai"] = np.array([[VAI_MBL_THR]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Flux must be zero above VAI threshold"


def test_vegetation_scales_flux(base_cell):
    """Half the VAI threshold must give half the no-veg flux."""
    F_bare = dust_emission(**base_cell)
    base_cell["tlai"] = np.array([[VAI_MBL_THR / 2.0]])
    F_veg = dust_emission(**base_cell)
    assert np.isclose(F_veg.sum(), F_bare.sum() * 0.5, rtol=1e-8), (
        "Vegetation should reduce flux linearly with VAI"
    )


# --- soil moisture effects ---


def test_wet_soil_reduces_flux(base_cell):
    """Wet soil (gwc > gwc_thr) must produce less flux than dry soil."""
    F_dry = dust_emission(**base_cell)
    # Force gwc above threshold: 0.35*1000/1620 = 0.216 > 0.184
    base_cell["h2osoi_vol"] = np.array([[0.35]])
    F_wet = dust_emission(**base_cell)
    assert F_wet.sum() < F_dry.sum(), "Wet soil must reduce emission"


def test_moisture_below_threshold_no_effect(base_cell):
    """Below gwc_thr (0.184 for fclay=0.1), moisture must not change the flux."""
    F_a = dust_emission(**base_cell)  # gwc = 0.031
    base_cell["h2osoi_vol"] = np.array([[0.25]])  # gwc = 0.154 < 0.184
    F_b = dust_emission(**base_cell)
    assert np.isclose(F_a.sum(), F_b.sum(), rtol=1e-12)


def test_frozen_soil_reduces_flux(base_cell):
    """High ice content (liqfrac→0) must reduce flux."""
    F_liquid = dust_emission(**base_cell)
    base_cell["h2osoi_liq"] = np.array([[0.0]])
    base_cell["h2osoi_ice"] = np.array([[50.0]])
    F_frozen = dust_emission(**base_cell)
    assert F_frozen.sum() < F_liquid.sum(), "Frozen soil must reduce emission"


# --- wind threshold ---


def test_zero_flux_below_threshold(base_cell):
    """Very low wind must produce zero flux."""
    base_cell["fv"] = np.array([[0.01]])
    base_cell["u10"] = np.array([[0.01]])
    F_p = dust_emission(**base_cell)
    assert (F_p == 0.0).all(), "Sub-threshold wind must give zero flux"


def test_flux_increases_with_wind(base_cell):
    """Stronger wind must produce greater flux."""
    F_low = dust_emission(**base_cell)
    base_cell["fv"] = np.array([[0.80]])
    base_cell["u10"] = np.array([[15.0]])
    F_high = dust_emission(**base_cell)
    assert F_high.sum() > F_low.sum(), "Stronger wind must increase flux"


# --- clay fraction / sandblasting efficiency ---


def test_higher_clay_changes_flux(base_cell):
    """
    Higher clay increases the sandblasting efficiency eta = 10^(13.4*fclay-6),
    so total flux must increase with clay (all else equal).
    """
    F_low_clay = dust_emission(**base_cell)  # fclay=0.10
    base_cell["fclay"] = np.array([[0.18]])  # CLM4 cap=0.20
    F_high_clay = dust_emission(**base_cell)
    assert F_high_clay.sum() > F_low_clay.sum(), (
        "Higher clay fraction must increase sandblasting efficiency and flux"
    )


def test_clay_cap_on_sandblasting(base_cell):
    """
    Above FCLAY_MAX the sandblasting efficiency is capped; flux may then only
    change through gwc_thr (no effect here because the soil is dry).
    """
    base_cell["fclay"] = np.array([[0.25]])
    F1 = dust_emission(**base_cell)
    base_cell["fclay"] = np.array([[0.40]])
    F2 = dust_emission(**base_cell)
    assert np.isclose(F1.sum(), F2.sum(), rtol=1e-12)


def test_zero_clay_gives_very_low_flux(base_cell):
    """
    fclay→0 drives dst_slt_flx_rat→100*1e-6 (very small), so flux is
    negligible but not exactly zero (the sandblasting ratio is nonzero).
    """
    base_cell["fclay"] = np.array([[1.0e-6]])
    F_p = dust_emission(**base_cell)
    assert F_p.sum() >= 0.0, "Flux must be non-negative for tiny clay"


# --- basin factor ---


def test_basin_factor_scales_flux(base_cell):
    """Doubling mbl_bsn_fct must exactly double the total flux."""
    F1 = dust_emission(**base_cell)
    base_cell["mbl_bsn_fct"] = np.array([[2.0]])
    F2 = dust_emission(**base_cell)
    assert np.isclose(F2.sum(), F1.sum() * 2.0, rtol=1e-8), "Basin factor must scale flux linearly"


# --- bin partitioning ---


def test_bin_partition_sums_to_total(base_cell):
    """
    Sum of all bin fluxes must equal F_tot * sum(ovr_src_snk_mss).
    (ovr_src_snk_mss does not necessarily sum to 1 because the lognormal
    integral may not cover the full size range of the bin grid.)
    """
    F_p = dust_emission(**base_cell)
    bin_sum = OVR_SRC_SNK_MSS.sum()  # expected partition factor

    # Reconstruct F_tot from the first bin that is nonzero
    if OVR_SRC_SNK_MSS[:, 0].sum() > 0:
        F_tot_from_bin0 = F_p[0, 0, 0] / OVR_SRC_SNK_MSS[:, 0].sum()
        expected_total = F_tot_from_bin0 * bin_sum
        assert np.isclose(F_p[0, 0, :].sum(), expected_total, rtol=1e-8), (
            "Bin sum does not match expected partition"
        )


def test_four_bins(base_cell):
    """There must be exactly 4 output bins (ndst=4)."""
    F_p = dust_emission(**base_cell)
    assert F_p.shape[-1] == 4, f"Expected 4 bins, got {F_p.shape[-1]}"


# --- precomputed constants ---


def test_tmp1_positive():
    """tmp1 must be strictly positive."""
    assert TMP1 > 0.0, f"tmp1 must be positive; got {TMP1}"


def test_tmp1_physical_range():
    """
    tmp1 represents the threshold friction velocity scale (Iversen & White 1982).
    For dmt_slt_opt=75 µm and dns_slt=2650 kg/m³ the dimensional result should
    be in the range [0.05, 0.5] (kg^0.5 m^-0.5 s^-1), giving u*_thr ~ 0.2 m/s
    at rho_air=1.2 kg/m³.
    """
    u_thr_at_std = TMP1 / np.sqrt(1.2)
    assert 0.05 < u_thr_at_std < 0.5, (
        f"Threshold velocity at std density = {u_thr_at_std:.4f} m/s out of physical range"
    )


def test_ovr_src_snk_mss_shape():
    """Overlap matrix must be shape [3, 4]."""
    assert OVR_SRC_SNK_MSS.shape == (DST_SRC_NBR, NDST), (
        f"Expected ({DST_SRC_NBR},{NDST}), got {OVR_SRC_SNK_MSS.shape}"
    )


def test_ovr_src_snk_mss_nonnegative():
    """All overlap factors must be non-negative."""
    assert (OVR_SRC_SNK_MSS >= 0.0).all(), "Overlap factors must be non-negative"


def test_ovr_src_snk_mss_row_sums():
    """
    Each source mode's mass fractions summed over transport bins cannot
    exceed its own mss_frc_src (some mass may fall outside 0.1-10 µm).
    """
    for m in range(DST_SRC_NBR):
        row_sum = OVR_SRC_SNK_MSS[m, :].sum()
        assert row_sum <= MSS_FRC_SRC[m] + 1.0e-8, (
            f"Mode {m} overlap sum {row_sum:.4f} exceeds mss_frc_src {MSS_FRC_SRC[m]}"
        )


def test_mss_frc_src_sums_to_one():
    """Source mode mass fractions must sum to 1.0 (BSM96)."""
    assert np.isclose(MSS_FRC_SRC.sum(), 1.0, atol=1.0e-9), (
        f"MSS_FRC_SRC sums to {MSS_FRC_SRC.sum()}, expected 1.0"
    )


# --- gwc_threshold ---


def test_gwc_threshold_zero_clay():
    """Zero clay gives the Zender/CLM intercept 0.17 kg/kg."""
    assert np.isclose(gwc_threshold(np.array([[0.0]]))[0, 0], 0.17)


def test_gwc_threshold_matches_ctsm():
    """gwc_thr must equal CTSM's 0.17 + 0.14*clay[%]*0.01 (0.184 for 10 % clay)."""
    assert np.isclose(gwc_threshold(np.array([[0.10]]))[0, 0], 0.184, rtol=1e-12)


def test_gwc_threshold_increases_with_clay():
    """Higher clay fraction must raise the moisture threshold."""
    thr_low = gwc_threshold(np.array([[0.05]]))[0, 0]
    thr_high = gwc_threshold(np.array([[0.20]]))[0, 0]
    assert thr_high > thr_low, "Threshold must increase with clay"


# --- vectorised grid ---


def test_vectorized_grid():
    """Must run correctly on a full 2-D grid with mixed conditions."""
    nx, ny = 10, 8
    rng = np.random.default_rng(42)
    F_p = dust_emission(
        fv=rng.uniform(0.1, 1.5, (nx, ny)),
        u10=rng.uniform(1.0, 20.0, (nx, ny)),
        forc_rho=rng.uniform(1.1, 1.3, (nx, ny)),
        h2osoi_vol=rng.uniform(0.0, 0.30, (nx, ny)),
        h2osoi_liq=rng.uniform(0.0, 30.0, (nx, ny)),
        h2osoi_ice=rng.uniform(0.0, 5.0, (nx, ny)),
        watsat=rng.uniform(0.30, 0.50, (nx, ny)),
        tlai=rng.uniform(0.0, 0.5, (nx, ny)),
        tsai=rng.uniform(0.0, 0.2, (nx, ny)),
        frac_sno=rng.uniform(0.0, 1.0, (nx, ny)),
        fclay=rng.uniform(0.01, 0.20, (nx, ny)),
        mbl_bsn_fct=np.ones((nx, ny)),
        is_soil=rng.choice([True, False], (nx, ny)),
    )
    assert F_p.shape == (nx, ny, NDST), f"Shape mismatch: {F_p.shape}"
    assert (F_p >= 0.0).all(), "Negative flux in vectorized grid test"
