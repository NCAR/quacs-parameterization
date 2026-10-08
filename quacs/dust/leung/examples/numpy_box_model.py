"""Standalone NumPy box-model demo for the leung dust emission scheme.

Run with: python examples/numpy_box_model.py
"""

import numpy as np

from quacs.dust.leung.leung import (
    BIT,
    CD0,
    CE,
    CK,
    DP_MED,
    DST_SRC_NBR,
    KAPPA_MAX,
    MSS_FRC_SRC,
    NDST,
    OVR_SRC_SNK_MSS,
    RHO_AIR0,
    U_FT0_DRY_REF,
    U_ST0,
    _log_law_factor,
    drag_partition,
    dust_emission,
    intermittency_factor,
    moisture_threshold_factor,
    sl2000_u_ft0_dry,
    turbulent_sigma_us,
)

if __name__ == "__main__":
    # ---- Representative semi-arid inputs -----------------------------------
    fv_val = 0.60  # [m/s]   friction velocity
    rho_val = 1.2  # [kg/m3] air density
    h2ovol_val = 0.05  # [m3/m3] volumetric soil water (dry)
    h2oliq_val = 10.0  # [kg/m2]
    h2oice_val = 0.0  # [kg/m2]
    watsat_val = 0.40  # [m3/m3]
    tlai_val = 0.0  # bare ground
    tsai_val = 0.0
    frac_sno_v = 0.0
    obu_val = -100.0  # [m]   unstable BL (convective)
    fclay_val = 0.10  # 10 % clay
    z0a_val = 1.0e-4  # [m]   smooth aeolian roughness
    Ar_val = 0.0  # no rocks
    Av_val = 0.0  # no plants
    Kc_val = 1.0  # unity (no upscaling correction)

    print("Leung et al. (2023/2024) Dust Emission Scheme — box model demo")
    print("=" * 68)
    print("Precomputed constants:")
    print(
        f"  u*ft0_dry at rho_air0 = {U_FT0_DRY_REF:.4f} m/s  "
        f"(Shao & Lu 2000, Dp = {DP_MED * 1e6:.0f} µm)"
    )
    print(f"  u*it  (Bit*u*ft0)     = {BIT * U_FT0_DRY_REF:.4f} m/s")
    print()
    print("Source-to-sink overlap matrix  ovr_src_snk_mss [3 modes × 4 bins]:")
    print("  Bin edges [µm]:  0.1-1.0  |  1.0-2.5  |  2.5-5.0  |  5.0-10.0")
    for m in range(DST_SRC_NBR):
        row = OVR_SRC_SNK_MSS[m, :]
        print(
            f"  Mode {m + 1} (mss_frc={MSS_FRC_SRC[m]:.3f}): "
            f"{row[0]:.5f}   {row[1]:.5f}   {row[2]:.5f}   {row[3]:.5f}"
        )
    print(
        "  Column sums (bin fractions): "
        + "   ".join(f"{OVR_SRC_SNK_MSS[:, n].sum():.5f}" for n in range(NDST))
    )
    print()

    F_p = dust_emission(
        fv=np.array([[fv_val]]),
        forc_rho=np.array([[rho_val]]),
        h2osoi_vol=np.array([[h2ovol_val]]),
        h2osoi_liq=np.array([[h2oliq_val]]),
        h2osoi_ice=np.array([[h2oice_val]]),
        watsat=np.array([[watsat_val]]),
        tlai=np.array([[tlai_val]]),
        tsai=np.array([[tsai_val]]),
        frac_sno=np.array([[frac_sno_v]]),
        obu=np.array([[obu_val]]),
        fclay=np.array([[fclay_val]]),
        z0a=np.array([[z0a_val]]),
        frac_rock=np.array([[Ar_val]]),
        frac_veg=np.array([[Av_val]]),
        Kc_map=np.array([[Kc_val]]),
        is_soil=np.array([[True]]),
    )

    # intermediate diagnostics
    fm, gwc_sfc, gwc_thr = moisture_threshold_factor(
        np.array([[h2ovol_val]]),
        np.array([[watsat_val]]),
        np.array([[fclay_val]]),
    )
    u_ft0_d = sl2000_u_ft0_dry(np.array([[rho_val]]))[0, 0]
    Feff_d, _, _ = drag_partition(
        np.array([[z0a_val]]),
        np.array([[0.0]]),
        np.array([[Ar_val]]),
        np.array([[Av_val]]),
    )
    u_s_d = fv_val * Feff_d[0, 0]
    u_it_d = BIT * u_ft0_d
    u_st_d = u_ft0_d * fm[0, 0] * np.sqrt(rho_val / RHO_AIR0)
    Cd_d = CD0 * np.exp(-CE * (u_st_d - U_ST0) / U_ST0)
    k_d = min(CK * (u_st_d - U_ST0) / U_ST0, KAPPA_MAX)
    sigma_d = turbulent_sigma_us(np.array([[u_s_d]]), np.array([[obu_val]]))[0, 0]
    eta_d = intermittency_factor(
        np.array([[u_s_d]]),
        np.array([[sigma_d]]),
        np.array([[u_it_d]]),
        np.array([[u_ft0_d * fm[0, 0]]]),
    )[0, 0]

    print(
        f"Cell diagnostics (fv={fv_val} m/s | rho={rho_val} kg/m3 | "
        f"fclay={fclay_val} | gwc_vol={h2ovol_val})"
    )
    print(
        f"  gwc_sfc = {gwc_sfc[0, 0]:.4f} kg/kg  |  "
        f"gwc_thr = {gwc_thr[0, 0]:.4f} kg/kg  |  fm = {fm[0, 0]:.4f}"
    )
    print(
        f"  u*ft0_dry = {u_ft0_d:.4f} m/s  |  u*it = {u_it_d:.4f} m/s  |  "
        f"u*ft_wet = {u_ft0_d * fm[0, 0]:.4f} m/s"
    )
    print(f"  Feff = {Feff_d[0, 0]:.4f}  |  u*s = {u_s_d:.4f} m/s  |  u*st = {u_st_d:.4f} m/s")
    print(
        f"  U(0.1 m) = {u_s_d * _log_law_factor():.3f} m/s  |  sigma_U(0.1 m) = {sigma_d:.4f} m/s"
    )
    print(f"  Cd = {Cd_d:.3e}  |  kappa = {k_d:.3f}  |  eta = {eta_d:.6f}")
    print()

    bin_labels = ["0.1-1 µm", "1-2.5 µm", "2.5-5 µm", "5-10 µm"]
    print("Per-bin flux [kg m-2 s-1]:")
    for n in range(NDST):
        print(f"  Bin {n + 1} ({bin_labels[n]}): {F_p[0, 0, n]:.4e}")
    print(f"  Total:             {F_p[0, 0, :].sum():.4e}")
    print()

    # flux-to-tendency time loop (mass mixing ratio)
    rho_air = rho_val
    dz = 50.0  # [m]  lowest model layer thickness
    dt = 300.0  # [s]  time step
    x = np.zeros(NDST)

    print(f"{'Step':>5}  {'F_total [kg/m2/s]':>20}  {'x_total [kg/kg]':>18}")
    print("-" * 52)
    for step in range(1, 11):
        tendency = F_p[0, 0, :] / (rho_air * dz)  # [kg kg-1 s-1]
        x += tendency * dt
        print(f"{step:>5}  {F_p[0, 0, :].sum():>20.4e}  {x.sum():>18.4e}")

    print()
    print("Run tests with: pytest quacs/dust/leung")
