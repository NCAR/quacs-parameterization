"""Standalone NumPy box-model demo for the kok dust emission scheme.

Run with: python examples/numpy_box_model.py
"""

import numpy as np

from quacs.dust.kok.kok import (
    C_KAPPA,
    C_TUNE_DEFAULT,
    CD0,
    CE,
    DP_DEFAULT,
    KAPPA_MAX,
    LAND,
    RHOP_DEFAULT,
    U_ST0,
    bare_fraction,
    dust_emission,
    erodibility_coefficient,
    fragmentation_exponent,
    moisture_factor,
    standardized_threshold,
    threshold_velocity_dry,
)

if __name__ == "__main__":
    # synthetic single-cell inputs
    ustar_val = 0.5  # m/s    friction velocity
    w_val = 0.02  # kg/kg  gravimetric soil moisture (dry)
    rho_a_val = 1.2  # kg/m^3 near-surface air density
    vai_val = 0.05  # [-]    sparse vegetation
    snowdp_val = 0.0  # m      no snow
    f_clay_val = 0.15  # [-]    typical desert clay fraction

    rho_a_1d = np.array([[rho_a_val]])
    f_clay_1d = np.array([[f_clay_val]])

    # --- threshold diagnostics ---
    u_ft0 = threshold_velocity_dry(Dp=DP_DEFAULT, rhop=RHOP_DEFAULT, rho_a=rho_a_1d)
    f_m = moisture_factor(np.array([[w_val]]), f_clay_1d)
    u_ft = u_ft0 * f_m
    u_st = standardized_threshold(u_ft, rho_a_1d)
    C_d = erodibility_coefficient(u_st)
    kappa = fragmentation_exponent(u_st)
    f_bare_val = bare_fraction(np.array([[vai_val]]))

    print("K14 dust emission scheme — box model demo")
    print("=" * 56)
    print("  Constants (Kok et al. 2014a):")
    print(f"    CD0     = {CD0:.2e}")
    print(f"    CE      = {CE}")
    print(f"    U_ST0   = {U_ST0} m/s")
    print(f"    C_kappa = {C_KAPPA}  (kappa_max = {KAPPA_MAX})")
    print(f"    C_tune  = {C_TUNE_DEFAULT}")
    print()
    print("  Threshold diagnostics:")
    print(f"    Dp      = {DP_DEFAULT * 1e6:.1f} µm")
    print(f"    u_ft0   = {u_ft0[0, 0]:.4f} m/s   (dry threshold, S&L00)")
    print(f"    f_m     = {f_m[0, 0]:.4f}          (Fécan moisture factor)")
    print(f"    u_ft    = {u_ft[0, 0]:.4f} m/s   (wet fluid threshold)")
    print(f"    u_st    = {u_st[0, 0]:.4f} m/s   (standardized threshold)")
    print()
    print("  Erodibility:")
    print(f"    C_d     = {C_d[0, 0]:.4e}       (soil erodibility coeff)")
    print(f"    kappa   = {kappa[0, 0]:.4f}         (fragmentation exponent)")
    print(f"    f_bare  = {f_bare_val[0, 0]:.4f}")
    print()

    # --- run a 10-step box model ---
    rho_air = 1.2  # kg/m^3
    dz = 50.0  # m      bottom layer thickness
    dt = 300.0  # s      model time step

    F_d = dust_emission(
        ustar=np.array([[ustar_val]]),
        w=np.array([[w_val]]),
        rho_a=rho_a_1d,
        vai=np.array([[vai_val]]),
        snowdp=np.array([[snowdp_val]]),
        f_clay=f_clay_1d,
        oro=np.array([[LAND]]),
    )

    x = 0.0  # tracer mixing ratio [kg/kg]
    print(
        f"  ustar={ustar_val} m/s | w={w_val} | f_clay={f_clay_val} | "
        f"vai={vai_val} | snowdp={snowdp_val}"
    )
    print()
    print(f"{'Step':>5}  {'F_d [kg/m2/s]':>18}  {'x [kg/kg]':>16}")
    print("-" * 45)
    for step in range(1, 11):
        tendency = F_d[0, 0] / (rho_air * dz)
        x += tendency * dt
        print(f"{step:>5}  {F_d[0, 0]:>18.4e}  {x:>16.4e}")

    print()
    print("Run tests with: pytest quacs/dust/kok")
