#!/usr/bin/env python3
"""
MUSICA box model: CESM2-SLH sea-salt bromine recycling in one marine cell.

Rates come from ``quacs.vsls.cesm2_slh.seasalt_rate_constants``; MICM integrates.
The run checks MICM against the analytic first-order solution and reports the
gas-phase bromine gain from sea-salt bromide (0.65 Br per reactant consumed).

Inputs below are illustrative marine-boundary-layer values, not observations.

Usage::

    python quacs/vsls/cesm2_slh/examples/musica_seasalt_box_model.py
"""

import numpy as np

from quacs.vsls.cesm2_slh import SEA_SALT_REACTIONS, run_box_model

MET = dict(
    temperature_k=288.0,
    pressure_pa=1.0e5,
    latitude_deg=40.0,
    calday=180.0,
    sad_seasalt_m2_m3=5.0e-5,  # 50 um2/cm3 (illustrative)
    ocean_ice_fraction=1.0,
)
INITIAL_PPTV = {"HOBr": 5.0, "BrONO2": 2.0, "BrNO2": 0.5, "Br2": 0.0, "BrCl": 0.0}
TIME_STEP_S = 600.0
N_STEPS = 144  # 24 h


def br_atoms(conc: dict) -> float:
    """Total gas-phase Br (mol/m3) over the tracked species."""
    return conc["HOBr"] + conc["BrONO2"] + conc["BrNO2"] + 2.0 * conc["Br2"] + conc["BrCl"]


def main():
    out = run_box_model(MET, INITIAL_PPTV, TIME_STEP_S, N_STEPS)
    rates, n_air = out["rates"], out["n_air"]
    c0 = {s: v[0] for s, v in out["concentrations"].items()}
    c = {s: v[-1] for s, v in out["concentrations"].items()}
    t_end = out["times"][-1]

    print("Sea-salt rate constants")
    for rxn in SEA_SALT_REACTIONS:
        k = rates[rxn.name]
        print(f"  {rxn.name}  {rxn.reactant:<7} k = {k:.3e} 1/s   tau = {1 / k / 3600:.1f} h")

    print(f"\nAfter {t_end / 3600:.0f} h (pptv): MICM vs analytic")
    consumed = 0.0
    for rxn in SEA_SALT_REACTIONS:
        exact = c0[rxn.reactant] * np.exp(-rates[rxn.name] * t_end)
        consumed += c0[rxn.reactant] - exact
        print(
            f"  {rxn.reactant:<7} {c[rxn.reactant] / n_air * 1e12:8.4f}"
            f"  {exact / n_air * 1e12:8.4f}"
        )
    for prod, yld in (("Br2", 0.65), ("BrCl", 0.35)):
        print(f"  {prod:<7} {c[prod] / n_air * 1e12:8.4f}  {yld * consumed / n_air * 1e12:8.4f}")

    gain = br_atoms(c) - br_atoms(c0)
    print(
        f"\nGas-phase Br gain from sea-salt bromide: {gain / n_air * 1e12:.4f} pptv"
        f"  (expected 0.65 x consumed = {0.65 * consumed / n_air * 1e12:.4f} pptv)"
    )


if __name__ == "__main__":
    main()
