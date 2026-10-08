"""Unit tests for the CAM wetdepa_v2 port (quacs.wetdep.aerosol_cam)."""

import numpy as np

from quacs.wetdep.aerosol_cam.wetdep_cam import (
    GRAVIT,
    cloud_volume_diag,
    column_burden,
    column_scavenging_coefficients,
    column_step_exponential,
    layer_scavenging_rates,
    wetdepa_v2,
)

# ---------------------------------------------------------------------------


def _single_cloud_column(
    nlev=6, cloud_k=1, cldfrac=0.4, precs_val=1.0e-7, evap_k=None, evap_frac=0.0
):
    """Helper: one stratiform cloud layer at cloud_k, optional evaporation layer."""
    shp = (1, nlev)
    pdel = np.full(shp, 5000.0)
    cldt = np.zeros(shp)
    cldt[0, cloud_k] = cldfrac
    cldc = np.zeros(shp)
    cwat = np.zeros(shp)
    cwat[0, cloud_k] = 2.0e-4
    precs = np.zeros(shp)
    precs[0, cloud_k] = precs_val
    evaps = np.zeros(shp)
    if evap_k is not None:
        flux = precs_val * 5000.0 / GRAVIT
        evaps[0, evap_k] = evap_frac * flux * GRAVIT / 5000.0
    zeros = np.zeros(shp)
    return dict(
        pdel=pdel,
        cldt=cldt,
        cldc=cldc,
        cwat=cwat,
        precs=precs,
        evaps=evaps,
        conicw=zeros,
        cmfdqr=zeros,
        evapc=zeros,
    )


def test_no_precip_no_scavenging():
    """Without precipitation production there is no wet removal."""
    shp = (2, 5)
    z = np.zeros(shp)
    r = wetdepa_v2(
        q=np.full(shp, 1e-8),
        pdel=np.full(shp, 5000.0),
        cldt=np.full(shp, 0.5),
        cldc=z,
        cwat=np.full(shp, 1e-4),
        precs=z,
        evaps=z,
        conicw=z,
        cmfdqr=z,
        evapc=z,
        deltat=1800.0,
        sol_fact=0.3,
    )
    assert np.all(r["scavt"] == 0.0)
    assert np.all(r["sfc_flux"] == 0.0)


def test_removal_nonpositive_resusp_nonnegative():
    col = _single_cloud_column(evap_k=4, evap_frac=0.5)
    r = wetdepa_v2(q=np.full((1, 6), 1e-8), deltat=1800.0, sol_fact=0.3, **col)
    assert np.all(r["ic_scavt"] <= 0.0) and np.all(r["bc_scavt"] <= 0.0)
    assert np.all(r["resusp"] >= 0.0)


def test_linear_in_tracer():
    """All terms are proportional to q: doubling q doubles the tendency."""
    col = _single_cloud_column(evap_k=4, evap_frac=0.3)
    q = np.linspace(1e-9, 5e-8, 6).reshape(1, 6)
    r1 = wetdepa_v2(q=q, deltat=1800.0, sol_fact=0.3, **col)
    r2 = wetdepa_v2(q=2.0 * q, deltat=1800.0, sol_fact=0.3, **col)
    assert np.allclose(r2["scavt"], 2.0 * r1["scavt"], rtol=1e-12, atol=0.0)


def test_column_mass_budget_explicit():
    """CAM form: -sum_k scavt*pdel/g equals the surface deposition flux."""
    col = _single_cloud_column(evap_k=4, evap_frac=0.4)
    q = np.full((1, 6), 3e-8)
    r = wetdepa_v2(q=q, deltat=1800.0, sol_fact=0.3, **col)
    lhs = -(r["scavt"] * col["pdel"] / GRAVIT).sum()
    assert np.isclose(lhs, r["sfc_flux"][0], rtol=1e-10)


def test_full_evaporation_resuspends_everything_from_above():
    """
    If all precipitation evaporates in layer 3, everything scavenged ABOVE
    layer 3 is released there.  Only what layer 3 itself scavenges (below-
    cloud, using the flux entering it) is carried further down, exactly as
    in CAM's recursion.
    """
    col = _single_cloud_column(evap_k=3, evap_frac=1.0)
    q = np.full((1, 6), 1e-8)
    r = wetdepa_v2(q=q, deltat=1800.0, sol_fact=0.3, **col)
    c, pdog = r["coef"], col["pdel"] / GRAVIT
    removed_above = ((c["k_strat"] + c["k_conv"])[0, :3] * q[0, :3] * pdog[0, :3]).sum()
    removed_at_3 = (c["k_strat"] + c["k_conv"])[0, 3] * q[0, 3] * pdog[0, 3]
    assert np.isclose(r["resusp"][0, 3] * pdog[0, 3], removed_above, rtol=1e-12)
    assert np.isclose(r["sfc_flux"][0], removed_at_3, rtol=1e-12)


def test_cldvst_below_single_cloud():
    """Below one precipitating cloud layer, the raining area equals its cloud fraction."""
    col = _single_cloud_column(cloud_k=1, cldfrac=0.4)
    cldvst, cldvcu = cloud_volume_diag(
        col["pdel"],
        col["cldt"],
        col["cldc"],
        col["precs"],
        col["evaps"],
        col["cmfdqr"],
        col["evapc"],
    )
    assert cldvst[0, 0] == 0.0 and cldvst[0, 1] == 0.0  # at/above the cloud
    assert np.allclose(cldvst[0, 2:], 0.4, rtol=1e-12)
    assert np.all(cldvcu == 0.0)


def test_below_cloud_rate_one_mm_per_hour():
    """
    Dimensional check of the Dana & Hales term: 1 mm/h of rain over the whole
    area (cldvst = 1), scavcoef = 0.1 mm-1, sol_fact = 1 -> k = 0.1 h-1.
    """
    precabs = 1.0 / 3600.0  # 1 mm/h = 1 kg m-2 h-1 -> kg m-2 s-1
    r = layer_scavenging_rates(
        cldt=0.0,
        cldc=0.0,
        cwat=0.0,
        precs=0.0,
        conicw=0.0,
        cmfdqr=0.0,
        precabs=precabs,
        precabc=0.0,
        cldvst=1.0,
        cldvcu=0.0,
        deltat=600.0,
        sol_fact=1.0,
        scavcoef=0.1,
    )
    assert np.isclose(r["k_total"], 0.1 / 3600.0, rtol=1e-12)


def test_in_cloud_rate_hand_calc():
    """k_ic = sol * clds * precs*dt/(cwat+precs*dt) / dt."""
    dt, sol, clds, cwat, precs = 1800.0, 0.3, 0.5, 2e-4, 1e-7
    r = layer_scavenging_rates(
        cldt=clds,
        cldc=0.0,
        cwat=cwat,
        precs=precs,
        conicw=0.0,
        cmfdqr=0.0,
        precabs=0.0,
        precabc=0.0,
        cldvst=0.0,
        cldvcu=0.0,
        deltat=dt,
        sol_fact=sol,
    )
    expected = sol * clds * (precs * dt / (cwat + precs * dt)) / dt
    assert np.isclose(r["k_total"], expected, rtol=1e-12)


def test_limiter_caps_removal():
    """Never remove more than the tracer in one step: k*dt <= 1."""
    r = layer_scavenging_rates(
        cldt=1.0,
        cldc=0.5,
        cwat=1e-6,
        precs=1e-3,
        conicw=1e-6,
        cmfdqr=1e-3,
        precabs=1.0,
        precabc=1.0,
        cldvst=1.0,
        cldvcu=1.0,
        deltat=1800.0,
        sol_fact=1.0,
    )
    assert r["k_total"] * 1800.0 <= 1.0


def test_exponential_step_conserves_mass():
    """Exact-in-time step: burden change == surface deposition * dt."""
    col = _single_cloud_column(evap_k=4, evap_frac=0.4)
    coef = column_scavenging_coefficients(deltat=1800.0, sol_fact=0.3, **col)
    q = np.full((1, 6), 3e-8)
    qn, s, flx = column_step_exponential(q, coef, 1800.0)
    dB = column_burden(q, col["pdel"]) - column_burden(qn, col["pdel"])
    assert np.isclose(dB[0], flx[0] * 1800.0, rtol=1e-10)


def test_exponential_matches_explicit_for_small_dt():
    """For k*dt -> 0 the exact step converges to CAM's explicit update."""
    col = _single_cloud_column(evap_k=4, evap_frac=0.4)
    dt = 1.0
    coef = column_scavenging_coefficients(deltat=dt, sol_fact=0.3, **col)
    q = np.full((1, 6), 3e-8)
    qn, _, _ = column_step_exponential(q, coef, dt)
    r = wetdepa_v2(q=q, deltat=dt, sol_fact=0.3, **col)
    assert np.allclose(qn, q + r["scavt"] * dt, rtol=1e-6, atol=1e-22)
