# -*- coding: utf-8 -*-
"""
musica_column_model.py
======================
MUSICA single-column model for CAM aerosol wet deposition
(wetdepa_v2 + clddiag, bulk / non-cloudborne branch).

The physics is in quacs/wetdep/aerosol_cam/wetdep_cam.py.

Each model level is one MICM grid cell.  Wet deposition couples the levels
through the precipitation: mass scavenged in upper layers falls with the
rain and is RESUSPENDED where the rain evaporates.  Per time step:

  1. Science module, top -> bottom sweep (column_step_exponential):
       k_k  = k_strat + k_conv                     [s-1]  first-order removal
       s_k  = resuspension source from the layers above, built from the
              mass actually removed above during this step   [kg kg-1 s-1]
  2. MICM, every level at once:
       dq_k/dt = -k_k q_k + s_k
       LOSS.dust_wetdep[k]  = k_k                  [s-1]
       EMIS.dust_resusp[k]  = s_k * rho_k / MW     [mol m-3 s-1]
  3. Checks:
       (a) MICM q_k(t+dt) vs the exact per-level solution
               q_new = s/k + (q - s/k) exp(-k dt)
       (b) column mass budget
               burden(t) - burden(t+dt) = surface wet deposition * dt
           with burden = sum_k q_k pdel_k / g          [kg m-2]

For comparison the same column is also run with CAM's own explicit update
(q <- q + scavt*dt from wetdepa_v2).  The two agree for k*dt << 1; the
difference shown at the end is CAM's O(k dt) time-discretization error.

LOSS:  state.set_user_defined_rate_parameters({"LOSS.dust_wetdep": [k_1, ..., k_n]})
EMIS:  state.set_user_defined_rate_parameters({"EMIS.dust_resusp": [e_1, ..., e_n]})

Units
-----
    q          [kg kg-1]      tracer mass mixing ratio (science module)
    C = q rho  [kg m-3]       mass concentration
    MICM       [mol m-3]      = C / MW   (MW is a placeholder)
    pdel       [Pa];  pdel/g  [kg m-2]
    precip     [kg m-2 s-1]   = mm s-1 of liquid water
    deposition [kg m-2 s-1];  cumulative [kg m-2]  (printed in mg m-2)

References:
    CAM wetdep.F90 (wetdepa_v2, clddiag): https://github.com/ESCOMP/CAM
    Rasch et al. (2000), Tellus B, 52, 1025-1056.
    Dana and Hales (1976), Atmos. Environ., 10, 45-50.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.wetdep.aerosol_cam import (
    GRAVIT,
    SCAVCOEF_DEFAULT,
    column_burden,
    column_scavenging_coefficients,
    column_step_exponential,
    wetdepa_v2,
)

R_DRY = 287.04  # [J kg-1 K-1]

# ---------------------------------------------------------------------------
# 1. Vertical grid (k = 0 is the model top, as in CAM)
# ---------------------------------------------------------------------------
NLEV = 20
p_int = np.linspace(10000.0, 100000.0, NLEV + 1)  # [Pa] interfaces, top -> bottom
pmid = 0.5 * (p_int[:-1] + p_int[1:])  # [Pa]
pdel = np.diff(p_int)  # [Pa] (4500 Pa each)

# Standard-atmosphere-like temperature, isothermal above the tropopause
T = np.maximum(288.15 * (pmid / 101325.0) ** 0.190263, 216.65)  # [K]
rho = pmid / (R_DRY * T)  # [kg m-3]
dz = pdel / (rho * GRAVIT)  # [m]

# ---------------------------------------------------------------------------
# 2. Clouds and precipitation (static over the run)
# ---------------------------------------------------------------------------
strat = (pmid >= 50000.0) & (pmid <= 80000.0)  # stratiform deck 500-800 hPa
conv = (pmid >= 40000.0) & (pmid <= 90000.0)  # light convection 400-900 hPa
evapz = pmid > 85000.0  # sub-cloud evaporation zone

cldc = np.where(conv, 0.05, 0.0)  # [-]
clds = np.where(strat, 0.60, 0.0)  # [-]
cldt = clds + cldc  # [-]
cwat = np.where(strat, 2.0e-4, 0.0)  # [kg kg-1]
precs = np.where(strat, 8.0e-8, 0.0)  # [kg kg-1 s-1]
conicw = np.where(conv, 5.0e-4, 0.0)  # [kg kg-1]
cmfdqr = np.where(conv, 2.0e-8, 0.0)  # [kg kg-1 s-1]

# Evaporation: 20 % of the precipitation entering each sub-cloud layer
EVAP_FRAC = 0.20
evaps = np.zeros(NLEV)
evapc = np.zeros(NLEV)
_flx_s = 0.0
_flx_c = 0.0
for k in range(NLEV):
    pdog = pdel[k] / GRAVIT
    if evapz[k]:
        evaps[k] = EVAP_FRAC * _flx_s / pdog
        evapc[k] = EVAP_FRAC * _flx_c / pdog
    _flx_s += (precs[k] - evaps[k]) * pdog
    _flx_c += (cmfdqr[k] - evapc[k]) * pdog
SFC_RAIN_MM_HR = (_flx_s + _flx_c) * 3600.0

# Tracer-specific factors (set by the calling aerosol module in CAM)
SOL_FACT = 0.3  # [-]  demo value, not a CAM default
SCAVCOEF = SCAVCOEF_DEFAULT  # [mm-1] 0.1, CAM value for non-modal aerosol

# ---------------------------------------------------------------------------
# 3. Initial dust profile: well mixed below 700 hPa, decaying above
# ---------------------------------------------------------------------------
Q_BL = 5.0e-8  # [kg kg-1] (~ 50-60 µg m-3)
q0 = np.where(pmid >= 70000.0, Q_BL, Q_BL * np.exp(-(70000.0 - pmid) / 10000.0))

# ---------------------------------------------------------------------------
# 4. Time stepping settings and scavenging coefficients
# ---------------------------------------------------------------------------
dt = 1800.0  # s  CAM physics step (rates are computed with it)
n_steps = 12  # 6 h


def as2d(a):
    """Add the ncol axis: the science module uses [ncol, nlev]."""
    return a[None, :]


met = dict(
    pdel=as2d(pdel),
    cldt=as2d(cldt),
    cldc=as2d(cldc),
    cwat=as2d(cwat),
    precs=as2d(precs),
    evaps=as2d(evaps),
    conicw=as2d(conicw),
    cmfdqr=as2d(cmfdqr),
    evapc=as2d(evapc),
)

coef = column_scavenging_coefficients(deltat=dt, sol_fact=SOL_FACT, scavcoef=SCAVCOEF, **met)
k_tot = (coef["k_strat"] + coef["k_conv"])[0]  # [s-1] per level

# ---------------------------------------------------------------------------
# 5. MICM mechanism: one tracer per cell, wet loss + resuspension source
# ---------------------------------------------------------------------------
DUST_MW = 1.0  # kg/mol, dimensional placeholder (MICM works in mol m-3)

dust = mc.Species(name="dust", molecular_weight_kg_mol=DUST_MW)
gas = mc.Phase(name="gas", species=[dust])

# Note: older MUSICA releases do not accept `products=` on FirstOrderLoss;
# the column mass budget is checked from burdens and the surface flux instead.
wetdep_rxn = mc.FirstOrderLoss(
    name="dust_wetdep",
    scaling_factor=1.0,
    reactants=[dust],
    gas_phase=gas,
)
resusp_rxn = mc.Emission(
    name="dust_resusp",
    scaling_factor=1.0,
    products=[dust],
    gas_phase=gas,
)
mechanism = mc.Mechanism(
    name="quacs_wetdep_cam_column",
    species=[dust],
    phases=[gas],
    reactions=[wetdep_rxn, resusp_rxn],
)

solver = musica.MICM(mechanism=mechanism, solver_type=musica.SolverType.rosenbrock_standard_order)
state = solver.create_state(number_of_grid_cells=NLEV)
state.set_conditions(list(T), list(pmid))
state.set_concentrations({"dust": list(q0 * rho / DUST_MW)})
state.set_user_defined_rate_parameters({"LOSS.dust_wetdep": list(k_tot)})

# ---------------------------------------------------------------------------
# 6. Print setup diagnostics
# ---------------------------------------------------------------------------
print("--- MUSICA wet-deposition column model (CAM wetdepa_v2 + clddiag) ---")
print(
    f"  {NLEV} levels, {p_int[0] / 100:.0f}-{p_int[-1] / 100:.0f} hPa, dt={dt:.0f} s, "
    f"n_steps={n_steps} ({n_steps * dt / 3600:.0f} h)"
)
print(
    f"  sol_fact={SOL_FACT}, scavcoef={SCAVCOEF} mm-1, evaporation "
    f"{EVAP_FRAC * 100:.0f}%/layer below 850 hPa"
)
print(f"  Surface precipitation: {SFC_RAIN_MM_HR:.3f} mm/h")
print()
print(
    f"{'k':>3} {'p[hPa]':>7} {'dz[m]':>6} {'cldt':>5} {'cldvst':>6} "
    f"{'precabs[mm/h]':>13} {'k_ic[s-1]':>10} {'k_bc[s-1]':>10} "
    f"{'k_tot[s-1]':>10} {'fracev':>6} {'q0[kg/kg]':>10}"
)
for k in range(NLEV):
    k_ic = coef["k_st_ic"][0, k] + coef["k_cv_ic"][0, k]
    k_bc = coef["k_st_bc"][0, k] + coef["k_cv_bc"][0, k]
    print(
        f"{k:>3} {pmid[k] / 100:>7.1f} {dz[k]:>6.0f} {cldt[k]:>5.2f} "
        f"{coef['cldvst'][0, k]:>6.3f} "
        f"{(coef['precabs'][0, k] + coef['precabc'][0, k]) * 3600:>13.4f} "
        f"{k_ic:>10.3e} {k_bc:>10.3e} {k_tot[k]:>10.3e} "
        f"{coef['fracev'][0, k]:>6.3f} {q0[k]:>10.3e}"
    )
print()

# ---------------------------------------------------------------------------
# 7. Time loop
# ---------------------------------------------------------------------------
RTOL_STEP = 1.0e-4  # MICM vs exact per-level solution (relative to column max)
RTOL_BUDGET = 1.0e-5  # burden change vs surface deposition (relative to burden);
# the exact sweep conserves to roundoff, so this measures
# MICM's truncation error summed over the column

q_micm = q0.copy()
q_cam = q0.copy()  # CAM explicit reference
cum_dep = 0.0  # [kg m-2]
cam_dep = 0.0  # [kg m-2]
B0 = column_burden(q0, pdel)

print(
    f"{'time(h)':>7}  {'burden[mg/m2]':>13}  {'wetdep cum[mg/m2]':>17}  "
    f"{'CAM expl. cum':>13}  {'step err':>9}  {'budget err':>10}  match"
)
print("-" * 86)
print(f"{0.0:>7.1f}  {B0 * 1e6:>13.5f}  {0.0:>17.5f}  {0.0:>13.5f}  {'-':>9}  {'-':>10}")

all_pass = True
max_step_err = 0.0
max_budget_err = 0.0
s_last = np.zeros(NLEV)

for step in range(n_steps):
    # (1) science module: exact top-down sweep -> resuspension sources s_k
    q_ref, s, flx = column_step_exponential(as2d(q_micm), coef, dt)
    q_ref, s, flx = q_ref[0], s[0], flx[0]
    s_last = s

    # (2) MICM: every level at once
    state.set_user_defined_rate_parameters(
        {"EMIS.dust_resusp": list(s * rho / DUST_MW)}
    )  # mol m-3 s-1
    B_old = column_burden(q_micm, pdel)
    solver.solve(state, dt)
    q_micm = np.array(state.get_concentrations()["dust"]) * DUST_MW / rho

    # (3) checks
    step_err = np.max(np.abs(q_micm - q_ref)) / np.max(np.abs(q_ref))
    budget_err = abs((B_old - column_burden(q_micm, pdel)) - flx * dt) / B_old
    max_step_err = max(max_step_err, step_err)
    max_budget_err = max(max_budget_err, budget_err)
    cum_dep += flx * dt

    # CAM explicit reference (q <- q + scavt*dt)
    rc = wetdepa_v2(q=as2d(q_cam), deltat=dt, sol_fact=SOL_FACT, scavcoef=SCAVCOEF, **met)
    q_cam = q_cam + rc["scavt"][0] * dt
    cam_dep += rc["sfc_flux"][0] * dt

    ok = (step_err <= RTOL_STEP) and (budget_err <= RTOL_BUDGET)
    all_pass = all_pass and ok
    print(
        f"{(step + 1) * dt / 3600:>7.1f}  {column_burden(q_micm, pdel) * 1e6:>13.5f}  "
        f"{cum_dep * 1e6:>17.5f}  {cam_dep * 1e6:>13.5f}  {step_err:>9.1e}  "
        f"{budget_err:>10.1e}  {'ok' if ok else 'FAIL'}"
    )

# ---------------------------------------------------------------------------
# 8. Final profile
# ---------------------------------------------------------------------------
print()
print(
    f"{'k':>3} {'p[hPa]':>7} {'q0':>11} {'q MICM':>11} {'q CAM expl':>11} "
    f"{'remain':>7} {'resusp[kg/kg/s]':>15}"
)
for k in range(NLEV):
    print(
        f"{k:>3} {pmid[k] / 100:>7.1f} {q0[k]:>11.4e} {q_micm[k]:>11.4e} "
        f"{q_cam[k]:>11.4e} {q_micm[k] / q0[k]:>7.3f} {s_last[k]:>15.4e}"
    )

B_end = column_burden(q_micm, pdel)
print()
print(
    f"  Column burden: {B0 * 1e6:.4f} -> {B_end * 1e6:.4f} mg m-2 "
    f"({(1 - B_end / B0) * 100:.1f} % removed in {n_steps * dt / 3600:.0f} h)"
)
print(
    f"  Cumulative wet deposition: MICM/exact {cum_dep * 1e6:.4f} mg m-2, "
    f"CAM explicit {cam_dep * 1e6:.4f} mg m-2 "
    f"({(cam_dep / cum_dep - 1) * 100:+.2f} %)"
)
print(f"  Mass closure: (B0 - B_end) - wetdep = {(B0 - B_end - cum_dep) * 1e6:.3e} mg m-2")
print(f"  max step error = {max_step_err:.2e},  max budget error = {max_budget_err:.2e}")
print()
if all_pass:
    print(
        f"All steps match the exact per-level solution (rtol={RTOL_STEP:.0e}) and "
        f"close the column mass budget (rtol={RTOL_BUDGET:.0e}). PASS"
    )
else:
    print("One or more steps failed the exact-solution or mass-budget check. FAIL")
