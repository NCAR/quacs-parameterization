#!/usr/bin/env python3
"""
MUSICA box model for Wesely (1989) gas dry deposition.

Computes the first-order loss rate k = Vd / H for each species with
``wesely_gas`` and integrates the concentrations with MICM.

Usage::

    python quacs/drydep/wesely/examples/musica_box_model.py
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.drydep.wesely import wesely_gas

BOX_HEIGHT_M = 1000.0  # m
SPECIES = ["O3", "SO2", "NO2", "HNO3", "NH3"]

MET = {
    "sfc_temp": 298.0,
    "air_temp": 296.0,
    "pressure_sfc": 101325.0,
    "pressure_10m": 100000.0,
    "wind_speed": 5.0,
    "spec_hum": 0.010,
    "solar_flux": 400.0,
    "month": 7,
    "lat": 40.0,
}
FRAC_LANDUSE = np.array([0.05, 0.30, 0.10, 0.25, 0.10, 0.05, 0.05, 0.00, 0.05, 0.05, 0.00])

if __name__ == "__main__":
    k_rates = {
        sp: wesely_gas(sp, frac_landuse=FRAC_LANDUSE, box_height_m=BOX_HEIGHT_M, **MET)["k"]
        for sp in SPECIES
    }
    for sp, k in k_rates.items():
        print(f"  {sp:<5}  v_d = {k * BOX_HEIGHT_M * 100:.4f} cm/s  k = {k:.4e} 1/s")

    species = [mc.Species(name=sp) for sp in SPECIES]
    gas = mc.Phase(name="gas", species=species)
    reactions = [
        mc.FirstOrderLoss(name=f"{s.name}_drydep", scaling_factor=1.0, reactants=[s], gas_phase=gas)
        for s in species
    ]
    mechanism = mc.Mechanism(
        name="wesely_drydep_box", species=species, phases=[gas], reactions=reactions
    )
    solver = musica.MICM(
        mechanism=mechanism, solver_type=musica.SolverType.rosenbrock_standard_order
    )

    state = solver.create_state(1)
    state.set_conditions(temperatures=[MET["air_temp"]], pressures=[MET["pressure_sfc"]])
    state.set_concentrations({sp: [1.0] for sp in SPECIES})
    state.set_user_defined_rate_parameters({f"LOSS.{sp}_drydep": [k_rates[sp]] for sp in SPECIES})

    time_step = 600.0  # s
    n_steps = 18  # 3 hours
    print(f"\n{'Time (h)':>8}  " + "  ".join(f"{sp:>10}" for sp in SPECIES))
    for i in range(n_steps + 1):
        if i > 0:
            solver.solve(state, time_step)
        concs = state.get_concentrations()
        row = "  ".join(f"{concs[sp][0]:10.6f}" for sp in SPECIES)
        print(f"{i * time_step / 3600:8.2f}  {row}")

    # MICM must agree with the exact solution C/C0 = exp(-k t).
    t_end = n_steps * time_step
    concs = state.get_concentrations()
    for sp in SPECIES:
        exact = np.exp(-k_rates[sp] * t_end)
        assert np.isclose(concs[sp][0], exact, rtol=1e-3), (sp, concs[sp][0], exact)
    print("\nMICM agrees with the exact exponential decay.")
