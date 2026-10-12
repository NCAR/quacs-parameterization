"""Standalone NumPy box-model demo for the zender dust emission scheme.

Run with: python examples/numpy_box_model.py
"""

import numpy as np

from quacs.dust.zender.zender import (
    DST_SRC_NBR,
    MSS_FRC_SRC,
    NDST,
    OVR_SRC_SNK_MSS,
    TMP1,
    dust_emission,
)

if __name__ == "__main__":
    # Synthetic single-cell inputs — representative semi-arid conditions
    fv_val = 0.50  # [m/s]    friction velocity
    u10_val = 8.0  # [m/s]    10-m wind
    rho_val = 1.2  # [kg/m3]  air density
    h2ovol_val = 0.05  # [m3/m3]  volumetric soil water (dry)
    h2oliq_val = 10.0  # [kg/m2]
    h2oice_val = 0.0  # [kg/m2]
    watsat_val = 0.40  # [m3/m3]
    tlai_val = 0.0  # bare ground
    tsai_val = 0.0
    frac_sno_val = 0.0
    fclay_val = 0.10  # 10 % clay
    mbl_bsn_val = 1.0

    print("DEAD Dust Emission (Zender et al. 2003a / CLM4) — box model demo")
    print("=" * 65)
    print("Precomputed constants:")
    print(f"  tmp1  = {TMP1:.6f}  [kg^0.5 m^-0.5 s^-1]")
    print(f"  u*_thr at rho=1.2 kg/m3 = {TMP1 / np.sqrt(rho_val):.4f} m/s")
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
        u10=np.array([[u10_val]]),
        forc_rho=np.array([[rho_val]]),
        h2osoi_vol=np.array([[h2ovol_val]]),
        h2osoi_liq=np.array([[h2oliq_val]]),
        h2osoi_ice=np.array([[h2oice_val]]),
        watsat=np.array([[watsat_val]]),
        tlai=np.array([[tlai_val]]),
        tsai=np.array([[tsai_val]]),
        frac_sno=np.array([[frac_sno_val]]),
        fclay=np.array([[fclay_val]]),
        mbl_bsn_fct=np.array([[mbl_bsn_val]]),
        is_soil=np.array([[True]]),
    )

    # flux-to-tendency time loop (QUACS doc v1, Section 7A)
    rho_air = rho_val
    dz = 50.0  # [m]  lowest model layer thickness
    dt = 300.0  # [s]  time step
    x = np.zeros(NDST)

    bin_labels = ["0.1-1 µm", "1-2.5 µm", "2.5-5 µm", "5-10 µm"]
    print(f"Cell: fv={fv_val} m/s | u10={u10_val} m/s | fclay={fclay_val} | gwet_vol={h2ovol_val}")
    print("Per-bin flux [kg m-2 s-1]:")
    for n in range(NDST):
        print(f"  Bin {n + 1} ({bin_labels[n]}): {F_p[0, 0, n]:.4e}")
    print(f"  Total:             {F_p[0, 0, :].sum():.4e}")
    print()
    print(f"{'Step':>5}  {'F_total [kg/m2/s]':>20}  {'x_total [kg/kg]':>18}")
    print("-" * 50)
    for step in range(1, 11):
        tendency = F_p[0, 0, :] / (rho_air * dz)  # [kg kg-1 s-1]
        x += tendency * dt
        print(f"{step:>5}  {F_p[0, 0, :].sum():>20.4e}  {x.sum():>18.4e}")

    print()
    print("Run tests with: pytest quacs/dust/zender")
