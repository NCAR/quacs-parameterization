#!/usr/bin/env python3
"""
MUSICA box model for Hg0 dry deposition.

Uses the offline GEOS-Chem dry deposition scheme (box_model.py) to
compute v_d (cm/s), converts it to a first-order loss rate k (s⁻¹) for a
1 m³ box, and integrates the concentration forward in time with MICM.

NOTE: this file is a work in progress. The mechanism definition below
contains a placeholder reaction ("A") that is not yet fully implemented.
"""

import numpy as np

import musica
import musica.mechanism_configuration as mc

from quacs.drydep.simple.box_model import DEFAULT_MET, XMW, compute_drydep_rate

# ── 1. Dry deposition rate ────────────────────────────────────────────────────
BOX_HEIGHT_M = 1.0  # m  (1 m³ box with 1 m² footprint)

dvel_cms, k_drydep = compute_drydep_rate(box_height_m=BOX_HEIGHT_M)

print(f"Dry deposition velocity:  {dvel_cms:.4f} cm s⁻¹")
print(f"First-order loss rate k:  {k_drydep:.4e} s⁻¹  (H = {BOX_HEIGHT_M} m)")
print(f"Hg0 e-folding lifetime:   {1/k_drydep:.1f} s  ({1/k_drydep/60:.1f} min)")


# ── 2. MUSICA mechanism ───────────────────────────────────────────────────────
Hg0_sp = mc.Species(name="Hg0", molecular_weight_kg_mol=XMW)
gas = mc.Phase(name="gas", species=[Hg0_sp])

dd_rxn = mc.FirstOrderLoss(
    name="Hg0_drydep",
    scaling_factor=1.0,
    reactants=[Hg0_sp],
    gas_phase=gas,
)

mechanism = mc.Mechanism(
    name="hg0_drydep_box",
    species=[Hg0_sp],
    phases=[gas],
    reactions=[dd_rxn],
)

solver = musica.MICM(
    mechanism=mechanism,
    solver_type=musica.SolverType.rosenbrock_standard_order,
)


# ── 3. Initial conditions ─────────────────────────────────────────────────────
# Typical atmospheric Hg0 ~ 1.5 ng m⁻³
Hg0_init_ng_m3 = 1.5
Hg0_init_mol_m3 = Hg0_init_ng_m3 * 1e-9 / (XMW * 1e3)  # ng/m³ → mol/m³

state = solver.create_state(1)
state.set_conditions(
    temperatures=[DEFAULT_MET["TC0"]],
    pressures=[DEFAULT_MET["PRESSU"]],
)
state.set_concentrations({"Hg0": [Hg0_init_mol_m3]})
state.set_user_defined_rate_parameters({"LOSS.Hg0_drydep": [k_drydep]})


# ── 4. Time integration ───────────────────────────────────────────────────────
# Run for 5 e-folding lifetimes so depletion is clearly visible.
time_step = 10.0
n_steps = int(5.0 / k_drydep / time_step)

times = np.zeros(n_steps + 1)
Hg0_concs = np.zeros(n_steps + 1)
Hg0_concs[0] = Hg0_init_mol_m3

for i in range(n_steps):
    solver.solve(state, time_step)
    Hg0_concs[i + 1] = state.get_concentrations()["Hg0"][0]
    times[i + 1] = (i + 1) * time_step

Hg0_ng_m3 = Hg0_concs * XMW * 1e12  # mol/m³ → ng/m³


# ── 5. Print results ──────────────────────────────────────────────────────────
print(f"\n{'Time (s)':>10}  {'Time (min)':>10}  {'Hg0 (ng/m3)':>14}  {'Fraction':>10}")
print("-" * 50)
for t, c, ng in zip(times, Hg0_concs, Hg0_ng_m3):
    print(f"{t:10.1f}  {t/60:10.2f}  {ng:14.6f}  {c/Hg0_init_mol_m3:10.6f}")
