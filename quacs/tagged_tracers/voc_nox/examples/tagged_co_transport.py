"""Run the tagged CO source-receptor case with MUSICA and 1-D transport.

Total CO has emissions in every grid cell. CO_tracer has emissions only in the
source cell. After each MUSICA step, an explicit 1-D advection-diffusion
operator moves both species. The script checks that CO_tracer gets to the
receptor cell, where its emission is zero.
"""

import musica
import numpy as np

from quacs.tagged_tracers.voc_nox import (
    build_1d_grid,
    build_co_tracer_mechanism,
    build_musica_rate_parameters,
    build_single_cell_tracer_emissions,
    build_total_co_emissions,
    derive_simple_met,
    run_source_receptor,
    source_receptor_proof_table,
)

NUM_GRID_CELLS = 30
SOURCE_CELL = 6
RECEPTOR_CELL = 20

DT = 10.0  # s
SIM_LENGTH = 3000.0  # s
DX = 1.0  # m
K_DIFF = 0.002  # m2 s-1
U_ADV = 0.006  # m s-1


def main():
    grid = build_1d_grid(NUM_GRID_CELLS, SOURCE_CELL, RECEPTOR_CELL)
    met = derive_simple_met(NUM_GRID_CELLS, temperature=298.0, pressure=101325.0)

    solver = musica.MICM(
        mechanism=build_co_tracer_mechanism(),
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(NUM_GRID_CELLS)
    state.set_conditions(temperatures=met["temperature"], pressures=met["pressure"])
    state.set_concentrations(
        {"CO": np.zeros(NUM_GRID_CELLS), "CO_tracer": np.zeros(NUM_GRID_CELLS)}
    )

    co_emis = build_total_co_emissions(grid)
    co_tracer_emis = build_single_cell_tracer_emissions(co_emis, SOURCE_CELL)
    rate_parameters = build_musica_rate_parameters(co_emis, co_tracer_emis)

    times = np.arange(0, SIM_LENGTH + DT, DT)
    co_hist, co_tracer_hist = run_source_receptor(
        solver, state, rate_parameters, times, dt=DT, K_diff=K_DIFF, u_adv=U_ADV, dx=DX
    )

    proof = source_receptor_proof_table(
        times, co_hist, co_tracer_hist, co_tracer_emis, SOURCE_CELL, RECEPTOR_CELL
    )
    print(proof.tail())

    assert co_tracer_emis[RECEPTOR_CELL] == 0.0
    assert proof["CO_tracer_at_receptor"].max() > 0.0
    print("PASSED: CO_tracer at the receptor is not zero, but its emission there is zero.")


if __name__ == "__main__":
    main()
