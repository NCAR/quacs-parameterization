# -*- coding: utf-8 -*-
"""
musica_box_model.py
===================
MUSICA box-model for CAM aerosol wet deposition (wetdepa_v2, bulk /
non-cloudborne branch).

The physics is in quacs/wetdep/aerosol_cam/wetdep_cam.py.

The box is ONE model layer inside a raining stratiform cloud, with a
prescribed precipitation flux falling into it from above.  The science
module returns first-order scavenging rate coefficients

    k = k_in-cloud + k_below-cloud                          [s-1]

which are handed to MICM as a FirstOrderLoss rate through the LOSS.*
user-defined parameter convention:

LOSS:  state.set_user_defined_rate_parameters({"LOSS.dust_wetdep": [k]})
EMIS:  state.set_user_defined_rate_parameters({"EMIS.dust_emis_total": [e]})

The mass removed by rain is diagnosed as  dust_total(0) - dust_total(t).
(Older MUSICA releases do not accept `products=` on FirstOrderLoss, so no
accumulator species is used.)

Units
-----
    k                   [s-1]          (from the science module)
    C (dust_total)      [kg m-3]       mass concentration used for output
    MICM concentration  [mol m-3]      = C / MW   (MW is a placeholder)
    precabs             [kg m-2 s-1]   = mm s-1 of liquid water
    scavcoef            [mm-1]         = m2 kg-1

Closed-form check: with constant k,
    C(t) = C0 * exp(-k t)             (MICM / exact)
CAM itself applies the tendency explicitly once per physics step:
    C_CAM(n dt) = C0 * (1 - k dt)^n
The two differ by O(k dt); both are printed.

References:
    CAM wetdep.F90 (wetdepa_v2, clddiag): https://github.com/ESCOMP/CAM
    Rasch et al. (2000), Tellus B, 52, 1025-1056.
    Dana and Hales (1976), Atmos. Environ., 10, 45-50.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.wetdep.aerosol_cam import (
    SCAVCOEF_DEFAULT,
    layer_scavenging_rates,
)

# ---------------------------------------------------------------------------
# 1. MICM mechanism: one tracer, one first-order loss with an accumulator
# ---------------------------------------------------------------------------
# Molecular weight is a dimensional placeholder (kg/mol).
# MICM works in mol m-3; all kg <-> mol conversions divide/multiply by DUST_MW.
DUST_MW = 1.0  # kg/mol

dust_total = mc.Species(name="dust_total", molecular_weight_kg_mol=DUST_MW)
gas = mc.Phase(name="gas", species=[dust_total])

# Note: older MUSICA releases do not accept `products=` on FirstOrderLoss,
# so the removed mass is diagnosed as C0 - C(t) instead of a product species.
wetdep_rxn = mc.FirstOrderLoss(
    name="dust_wetdep",
    scaling_factor=1.0,
    reactants=[dust_total],
    gas_phase=gas,
)

mechanism = mc.Mechanism(
    name="quacs_wetdep_cam_box",
    species=[dust_total],
    phases=[gas],
    reactions=[wetdep_rxn],
)

# ---------------------------------------------------------------------------
# 2. Layer meteorology (single layer inside a raining stratiform cloud)
# ---------------------------------------------------------------------------
T_LAYER = 270.0  # K
P_LAYER = 70000.0  # Pa  (~700 hPa)

cldt = 0.5  # [-]           total cloud fraction
cldc = 0.0  # [-]           convective cloud fraction (none here)
cwat = 2.0e-4  # [kg kg-1]     stratiform cloud water (0.2 g/kg)
precs = 1.0e-7  # [kg kg-1 s-1] stratiform precip production in this layer
conicw = 0.0  # [kg kg-1]     convective cloud water
cmfdqr = 0.0  # [kg kg-1 s-1] convective precip production

# Precipitation falling INTO the layer from above (a box has no column):
RAIN_MM_PER_HR = 1.0
precabs = RAIN_MM_PER_HR / 3600.0  # [kg m-2 s-1]  (1 mm h-1 = 1 kg m-2 h-1)
precabc = 0.0  # [kg m-2 s-1]
cldvst = 0.5  # [-]  raining area of the precipitation from above
cldvcu = 0.0  # [-]

# Tracer-specific factors (set by the calling aerosol module in CAM)
SOL_FACT = 0.3  # [-]  solubility factor — demo value, not a CAM default
SCAVCOEF = SCAVCOEF_DEFAULT  # [mm-1] 0.1, CAM value for non-modal aerosol

# ---------------------------------------------------------------------------
# 3. Integration settings (CAM physics time step)
# ---------------------------------------------------------------------------
dt = 1800.0  # s  (rates are computed with the physics step, as in CAM)
n_steps = 12  # 12 * 1800 s = 6 h
C0 = 1.0e-7  # [kg m-3]  initial dust (100 µg m-3 plume aloft)

# ---------------------------------------------------------------------------
# 4. Scavenging rate coefficients from the science module
# ---------------------------------------------------------------------------
r = layer_scavenging_rates(
    cldt=cldt,
    cldc=cldc,
    cwat=cwat,
    precs=precs,
    conicw=conicw,
    cmfdqr=cmfdqr,
    precabs=precabs,
    precabc=precabc,
    cldvst=cldvst,
    cldvcu=cldvcu,
    deltat=dt,
    sol_fact=SOL_FACT,
    scavcoef=SCAVCOEF,
)
k_total = float(r["k_total"])  # [s-1]

print("--- MUSICA wet-deposition box model (CAM wetdepa_v2, bulk aerosol) ---")
print(
    f"  Layer: T={T_LAYER:.0f} K, p={P_LAYER / 100:.0f} hPa, cldt={cldt}, "
    f"cwat={cwat:.1e} kg/kg, precs={precs:.1e} kg/kg/s"
)
print(
    f"  Rain from above: {RAIN_MM_PER_HR:.1f} mm/h  (precabs={precabs:.4e} kg m-2 s-1), "
    f"raining area cldvst={cldvst}"
)
print(f"  sol_fact={SOL_FACT},  scavcoef={SCAVCOEF} mm-1,  dt={dt:.0f} s")
print()
print("  Scavenging diagnostics:")
print(f"    fracp (strat, in-cloud)  = {float(r['fracp_s']):.4f}   [precs*dt/(cwat+precs*dt)]")
print(f"    odds  (strat, below)     = {float(r['odds_s']):.4f}   [precabs/cldvst*scavcoef*dt]")
print(f"    k in-cloud               = {float(r['k_st_ic']):.4e} s-1")
print(f"    k below-cloud            = {float(r['k_st_bc']):.4e} s-1")
print(
    f"    k total (after limiter)  = {k_total:.4e} s-1   "
    f"(e-folding time {1.0 / k_total / 3600.0:.2f} h,  k*dt = {k_total * dt:.4f})"
)
print()

# ---------------------------------------------------------------------------
# 5. Build solver and initial state
# ---------------------------------------------------------------------------
solver = musica.MICM(
    mechanism=mechanism,
    solver_type=musica.SolverType.rosenbrock_standard_order,
)
state = solver.create_state(number_of_grid_cells=1)
state.set_conditions([T_LAYER], [P_LAYER])
state.set_concentrations({"dust_total": [C0 / DUST_MW]})
state.set_user_defined_rate_parameters({"LOSS.dust_wetdep": [k_total]})

# ---------------------------------------------------------------------------
# 6. Time loop
# ---------------------------------------------------------------------------
times = np.zeros(n_steps + 1)
conc = np.zeros(n_steps + 1)
conc[0] = C0  # [kg m-3]
dep = np.zeros(n_steps + 1)  # [kg m-3] removed by rain

for step in range(n_steps):
    solver.solve(state, dt)
    c = state.get_concentrations()
    conc[step + 1] = c["dust_total"][0] * DUST_MW  # mol m-3 -> kg m-3
    dep[step + 1] = C0 - conc[step + 1]  # removed by rain [kg m-3]
    times[step + 1] = (step + 1) * dt

# ---------------------------------------------------------------------------
# 7. Print evolution + checks
# ---------------------------------------------------------------------------
RTOL_EXACT = 1.0e-4  # MICM vs exp(-k t): Rosenbrock truncation error at k*dt ~ 0.09
# (an emission-only box is exact; a loss term is not)
print(
    f"{'time (h)':>8}  {'dust [kg/m3]':>14}  {'exact [kg/m3]':>14}  "
    f"{'CAM explicit':>14}  {'deposited':>12}  {'rel err':>9}  match"
)
print("-" * 92)

all_pass = True
max_rel_err = 0.0
for n in range(n_steps + 1):
    exact = C0 * np.exp(-k_total * times[n])
    cam_expl = C0 * (1.0 - k_total * dt) ** n
    rel_err = abs(conc[n] - exact) / exact
    max_rel_err = max(max_rel_err, rel_err)

    ok = rel_err <= RTOL_EXACT
    all_pass = all_pass and ok
    print(
        f"{times[n] / 3600:>8.1f}  {conc[n]:>14.6e}  {exact:>14.6e}  "
        f"{cam_expl:>14.6e}  {dep[n]:>12.4e}  {rel_err:>9.1e}  "
        f"{'ok' if ok else 'FAIL'}"
    )

print()
print(f"  max |MICM - exact| / exact = {max_rel_err:.2e}")
print(
    f"  remaining after {times[-1] / 3600:.0f} h: MICM {conc[-1] / C0:.4f}, "
    f"CAM explicit {(1.0 - k_total * dt) ** n_steps:.4f}  (difference is the O(k dt) "
    f"time-discretization error of CAM's explicit update)"
)
print()
if all_pass:
    print(f"All steps match exp(-k t) (rtol={RTOL_EXACT:.0e}). PASS")
else:
    print("One or more steps deviated from exp(-k t). FAIL")
