# -*- coding: utf-8 -*-
"""
musica_box_model.py
===================
MUSICA box-model for the Kok et al. (2014) / K14 dust emission scheme.

The physics lives in quacs/dust/kok/kok.py (untouched science).

This file builds the MICM mechanism, creates the solver, and time-steps.
Each step it asks the science module for F_d,
divides by dz to get a volumetric emission rate,
and sets it on the MICM state via the EMIS.* user-defined parameter convention.

LOSS:  state.set_user_defined_rate_parameters({"LOSS.Hg0_drydep": [k]})
EMIS:  state.set_user_defined_rate_parameters({"EMIS.dust_emis_total": [e]})

Key difference from Ginoux / Zender box models
-----------------------------------------------
The K14 scheme (quacs.dust.kok.dust_emission) returns a single integrated
total flux F_d [x, y] in [kg m-2 s-1], not a per-bin array. The brittle-
fragmentation size distribution (Kok 2011) that partitions F_d into Aitken /
accumulation / coarse modes is applied at a later stage (CAM6 microphysics).
Therefore this box model uses ONE emitted species ("dust_total") and one
EMIS reaction.

The static and synthetic-dynamic inputs below mirror the standalone demo in
examples/numpy_box_model.py: same numbers, same threshold computation, but solved by
the Rosenbrock integrator instead of an explicit Forward Euler step.

References:
    Kok, J. F., et al. (2014). An improved dust emission model – Part 1.
    Atmos. Chem. Phys., 14, 13023–13041.
    https://doi.org/10.5194/acp-14-13023-2014

    Kok, J. F., et al. (2014). An improved dust emission model – Part 2.
    Atmos. Chem. Phys., 14, 13043–13061.
    https://doi.org/10.5194/acp-14-13043-2014

    Shao, Y. and Lu, H. (2000). A simple expression for wind erosion
    threshold friction velocity. J. Geophys. Res., 105, 22437–22443.
    https://doi.org/10.1029/2000JD900304

    Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization
    of the increase of the aeolian erosion threshold wind friction velocity.
    Ann. Geophys., 17, 149–157.
    https://doi.org/10.1007/s00585-999-0149-7
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.kok.kok import (
    C_KAPPA,
    C_TUNE_DEFAULT,
    CD0,
    CE,
    DP_DEFAULT,
    FECAN_A_DEFAULT,
    KAPPA_MAX,
    LAND,
    RHOP_DEFAULT,
    U_ST0,
    VAI_THR,
    bare_fraction,
    dust_emission,
    erodibility_coefficient,
    fragmentation_exponent,
    moisture_factor,
    standardized_threshold,
    threshold_velocity_dry,
)

# ---------------------------------------------------------------------------
# 1. MICM mechanism: single emitted species for total dust flux
#    (K14 outputs [x, y] total flux; bin partitioning is done downstream
#    by the brittle-fragmentation size distribution in CAM6 microphysics)
# ---------------------------------------------------------------------------
# Molecular weight is a dimensional placeholder (kg/mol).
# MICM works in mol m-3 (concentrations) and mol m-3 s-1 (EMIS rates).
# All kg <-> mol conversions below divide/multiply by DUST_MW explicitly,
# so results stay correct if this placeholder is ever changed.
DUST_MW = 1.0  # kg/mol

dust_total = mc.Species(name="dust_total", molecular_weight_kg_mol=DUST_MW)
gas = mc.Phase(name="gas", species=[dust_total])

dust_emis_rxn = mc.Emission(
    name="dust_emis_total",
    scaling_factor=1.0,
    products=[dust_total],
    gas_phase=gas,
)

mechanism = mc.Mechanism(
    name="quacs_dust_kok2014",
    species=[dust_total],
    phases=[gas],
    reactions=[dust_emis_rxn],
)

# ---------------------------------------------------------------------------
# 2. Static / synthetic inputs (single cell [1, 1])
#    Mirrors the demo in examples/numpy_box_model.py.
# ---------------------------------------------------------------------------
# Static land-surface properties
f_clay = np.array([[0.15]])  # clay mass fraction [-], 15 % (typical desert)
oro = np.array([[LAND]])  # land mask (1.0 = land)

# Vegetation / snow (static for this demo; dynamic in full model)
vai = np.array([[0.05]])  # vegetation area index [-], sparse veg
snowdp = np.array([[0.0]])  # snow depth [m], no snow

# ---------------------------------------------------------------------------
# 3. Pre-compute and print threshold diagnostics (mirrors examples/numpy_box_model.py)
# ---------------------------------------------------------------------------
rho_a_demo = np.array([[1.2]])
w_demo = np.array([[0.02]])

u_ft0_demo = threshold_velocity_dry(Dp=DP_DEFAULT, rhop=RHOP_DEFAULT, rho_a=rho_a_demo)
f_m_demo = moisture_factor(w_demo, f_clay)
u_ft_demo = u_ft0_demo * f_m_demo
u_st_demo = standardized_threshold(u_ft_demo, rho_a_demo)
C_d_demo = erodibility_coefficient(u_st_demo)
kappa_demo = fragmentation_exponent(u_st_demo)
f_bare_demo = bare_fraction(vai)

print("--- MUSICA dust-emission box model (Kok et al. 2014 / K14) ---")
print("  Constants (Kok et al. 2014a, Table 1):")
print(f"    CD0      = {CD0:.2e}")
print(f"    CE       = {CE}")
print(f"    U_ST0    = {U_ST0} m/s")
print(f"    C_kappa  = {C_KAPPA}  (kappa_max = {KAPPA_MAX})")
print(f"    C_tune   = {C_TUNE_DEFAULT}")
print()
print("  Threshold diagnostics (at demo forcing):")
print(f"    Dp       = {DP_DEFAULT * 1e6:.1f} µm")
print(f"    u_ft0    = {u_ft0_demo[0, 0]:.4f} m/s   (dry threshold, Shao & Lu 2000)")
print(f"    f_m      = {f_m_demo[0, 0]:.4f}          (Fécan moisture factor)")
print(f"    u_ft     = {u_ft_demo[0, 0]:.4f} m/s   (wet fluid threshold)")
print(f"    u_st     = {u_st_demo[0, 0]:.4f} m/s   (standardized threshold)")
print(f"    C_d      = {C_d_demo[0, 0]:.4e}       (soil erodibility coefficient)")
print(f"    kappa    = {kappa_demo[0, 0]:.4f}         (fragmentation exponent)")
print(f"    f_bare   = {f_bare_demo[0, 0]:.4f}")
print()


# ---------------------------------------------------------------------------
# 4. Helper: surface flux -> MICM emission rate
#    F_d [kg m-2 s-1] / dz [m]      = e_mass [kg m-3 s-1]
#    e_mass / DUST_MW [kg mol-1]    = e      [mol m-3 s-1]  (MICM source, dC/dt = +e)
# ---------------------------------------------------------------------------
def compute_emission_rate(ustar, w, rho_a, dz):
    """
    Call quacs.dust.kok.dust_emission for the current time step and convert
    the total surface flux [kg m-2 s-1] to a volumetric molar emission rate
    [mol m-3 s-1] for the MICM EMIS.dust_emis_total parameter.

    Parameters
    ----------
    ustar : np.ndarray [x, y]   friction velocity [m/s]
    w     : np.ndarray [x, y]   gravimetric soil moisture [kg/kg]
    rho_a : np.ndarray [x, y]   near-surface air density [kg/m³]
    dz    : float               lowest model layer thickness [m]

    Returns
    -------
    dict  {"EMIS.dust_emis_total": [e]}  with e in [mol m-3 s-1]
    """
    F_d = dust_emission(
        ustar=ustar,
        w=w,
        rho_a=rho_a,
        vai=vai,
        snowdp=snowdp,
        f_clay=f_clay,
        oro=oro,
        Dp=DP_DEFAULT,
        rhop=RHOP_DEFAULT,
        a=FECAN_A_DEFAULT,
        vai_thr=VAI_THR,
        C_tune=C_TUNE_DEFAULT,
    )  # [x, y]
    e_mass = F_d[0, 0] / dz  # [kg m-3 s-1]
    return {"EMIS.dust_emis_total": [e_mass / DUST_MW]}  # [mol m-3 s-1]


# ---------------------------------------------------------------------------
# 5. Build solver and initial state
# ---------------------------------------------------------------------------
solver = musica.MICM(
    mechanism=mechanism,
    solver_type=musica.SolverType.rosenbrock_standard_order,
)

state = solver.create_state(number_of_grid_cells=1)

state.set_conditions([300.0], [101325.0])  # T [K], p [Pa]

state.set_concentrations({"dust_total": [0.0]})

# ---------------------------------------------------------------------------
# 6. Synthetic dynamic forcing (constant over all steps; mirrors examples/numpy_box_model.py)
# ---------------------------------------------------------------------------
dt = 120.0  # s, model time step
dz = 50.0  # m, lowest-layer thickness (from MPAS-A / CLM5)
n_steps = 30  # 30 * 120 s = 1 hour

ustar_val = np.array([[0.5]])  # m/s, friction velocity (synthetic, constant)
w_val = np.array([[0.02]])  # kg/kg, gravimetric soil moisture (dry)
rho_a_val = np.array([[1.2]])  # kg/m³, near-surface air density

print("  Dynamic forcing (constant):")
print(
    f"    ustar = {ustar_val[0, 0]:.2f} m/s,  w = {w_val[0, 0]:.3f} kg/kg,  "
    f"rho_a = {rho_a_val[0, 0]:.1f} kg/m³"
)
print(f"    f_clay = {f_clay[0, 0]:.2f},  vai = {vai[0, 0]:.2f},  snowdp = {snowdp[0, 0]:.2f} m")
print(f"  Integration: dz={dz:.1f} m,  dt={dt:.0f} s,  n_steps={n_steps}\n")

# ---------------------------------------------------------------------------
# 7. Pre-compute emission rate (constant forcing → constant rate)
# ---------------------------------------------------------------------------
rates = compute_emission_rate(ustar_val, w_val, rho_a_val, dz)
state.set_user_defined_rate_parameters(rates)

e_total = rates["EMIS.dust_emis_total"][0]  # [mol m-3 s-1]
e_mass = e_total * DUST_MW  # [kg m-3 s-1]
print(f"  Surface flux F_d     : {e_mass * dz:.4e} kg m-2 s-1")
print(
    f"  Total emission rate  : {e_mass:.4e} kg m-3 s-1  "
    f"(= {e_total:.4e} mol m-3 s-1 with MW = {DUST_MW} kg/mol)\n"
)

# ---------------------------------------------------------------------------
# 8. Time loop: integrate with Rosenbrock solver
# ---------------------------------------------------------------------------
times = np.zeros(n_steps + 1)
concs = np.zeros(n_steps + 1)  # [kg m-3], single species (converted from mol m-3)

for step in range(n_steps):
    solver.solve(state, dt)

    state_conc = state.get_concentrations()
    concs[step + 1] = state_conc["dust_total"][0] * DUST_MW  # mol m-3 -> kg m-3
    times[step + 1] = (step + 1) * dt

# ---------------------------------------------------------------------------
# 9. Print evolution + closed-form check at every step
#    For a pure emission system with constant rate, C(t) = e * t exactly.
# ---------------------------------------------------------------------------
print(f"{'time (min)':>10}  {'dust_total [kg/m3]':>22}  {'expected [kg/m3]':>22}  match")
print("-" * 70)

all_pass = True

for step in range(n_steps + 1):
    actual = concs[step]
    expected = e_mass * times[step]  # [kg m-3]

    ok = np.isclose(actual, expected, rtol=1e-10)
    all_pass = all_pass and ok

    print(
        f"{times[step] / 60:>10.2f}  {actual:>22.6e}  {expected:>22.6e}  {'ok' if ok else 'FAIL'}"
    )

print()
if all_pass:
    print("All steps match the closed-form solution (rtol=1e-10). PASS")
else:
    print("One or more steps deviated from the closed-form solution. FAIL")
