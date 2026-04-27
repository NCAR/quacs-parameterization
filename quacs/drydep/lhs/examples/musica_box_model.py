#!/usr/bin/env python3
"""
MUSICA box model for multi-species dry deposition.

Computes v_d (cm/s) for Hg0, SO2, and O3 and integrates concentrations
forward in time with MICM over a tropical mixed-forest grid cell.

Usage::

    python -m quacs.drydep.lhs.examples.musica_box_model
"""

import numpy as np
import musica
import musica.mechanism_configuration as mc

from quacs.drydep.lhs import olson_land_cover
from quacs.drydep.lhs.box_model import compute_drydep_rate
from quacs.drydep.lhs.species import hg0, o3, so2

species = [hg0, so2, o3]
BOX_HEIGHT_M = 1.0  # m

met = dict(
    TC0=302.0,
    CFRAC=0.8,
    RADIAT=606.0,
    AZO=1.5,
    USTAR=0.15,
    PRESSU=1e5,
    SUNCOS_MID=0.8,
    ALBD=0.1,
    BXHEIGHT=130.0,
    U10M=0.1,
    V10M=-0.05,
    AIRDEN=1.15,
    HFLUX=16.8,
)

land_cover = [
    olson_land_cover.patch_from_olson(0, fraction=0.312, lai=0.10),
    olson_land_cover.patch_from_olson(29, fraction=0.200, lai=0.71),
    olson_land_cover.patch_from_olson(33, fraction=0.175, lai=0.61),
    olson_land_cover.patch_from_olson(43, fraction=0.313, lai=0.74),
]

if __name__ == "__main__":
    k_rates = compute_drydep_rate(
        species,
        box_height_m=BOX_HEIGHT_M,
        met=met,
        land_cover=land_cover,
        coefficients=olson_land_cover.coefficients,
    )

    print("Dry deposition velocities and loss rates")
    print("=" * 60)
    for sp in species:
        k = k_rates[sp.name][0]
        dvel_cms = k * BOX_HEIGHT_M / 0.01
        lifetime_s = 1.0 / k
        print(
            f"  {sp.name:<5}  v_d = {dvel_cms:.4f} cm s⁻¹  "
            f"k = {k:.4e} s⁻¹  τ = {lifetime_s:.1f} s ({lifetime_s/3600:.2f} h)"
        )
    print()

    gas = mc.Phase(name="gas", species=species)
    reactions = [
        mc.FirstOrderLoss(name=f"{sp.name}_drydep", scaling_factor=1.0,
                          reactants=[sp], gas_phase=gas)
        for sp in species
    ]
    mechanism = mc.Mechanism(name="drydep_multispecies_box", species=species,
                             phases=[gas], reactions=reactions)
    solver = musica.MICM(mechanism=mechanism,
                         solver_type=musica.SolverType.rosenbrock_standard_order)

    initial_ng_m3 = {"Hg0": 1.5, "SO2": 1000.0, "O3": 120e6}
    initial_mol_m3 = {
        sp.name: initial_ng_m3[sp.name] * 1e-9 / (sp.molecular_weight_kg_mol * 1e3)
        for sp in species
    }

    state = solver.create_state(1)
    state.set_conditions(temperatures=[met["TC0"]], pressures=[met["PRESSU"]])
    state.set_concentrations({sp.name: [initial_mol_m3[sp.name]] for sp in species})
    state.set_user_defined_rate_parameters(
        {f"LOSS.{sp.name}_drydep": [k_rates[sp.name][0]] for sp in species}
    )

    k_max = max(k_rates[sp.name][0] for sp in species)
    time_step = 10.0
    n_steps = int(3.0 / k_max / time_step)

    times = np.zeros(n_steps + 1)
    concs = {sp.name: np.zeros(n_steps + 1) for sp in species}
    for sp in species:
        concs[sp.name][0] = initial_mol_m3[sp.name]

    for i in range(n_steps):
        solver.solve(state, time_step)
        cell_concs = state.get_concentrations()
        for sp in species:
            concs[sp.name][i + 1] = cell_concs[sp.name][0]
        times[i + 1] = (i + 1) * time_step

    header = f"{'Time (s)':>10}  {'Time (min)':>10}  " + "  ".join(
        f"{sp.name:>14}" for sp in species
    )
    print(header)
    print("-" * len(header))
    for step in range(0, n_steps + 1, max(1, n_steps // 20)):
        t = times[step]
        row = f"{t:10.1f}  {t/60:10.2f}  "
        row += "  ".join(
            f"{concs[sp.name][step] * sp.molecular_weight_kg_mol * 1e12:14.6f}"
            for sp in species
        )
        print(row)
