"""Smoke test for the tagged tracer prototype.

This test only checks that the functions run. The scientist has not written unit tests yet.
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
)


def test_run_source_receptor_runs():
    n = 10
    grid = build_1d_grid(n, source_cell=2, receptor_cell=7)
    met = derive_simple_met(n, temperature=298.0, pressure=101325.0)
    solver = musica.MICM(
        mechanism=build_co_tracer_mechanism(),
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(n)
    state.set_conditions(temperatures=met["temperature"], pressures=met["pressure"])
    state.set_concentrations({"CO": np.zeros(n), "CO_tracer": np.zeros(n)})
    co_emis = build_total_co_emissions(grid)
    co_tracer_emis = build_single_cell_tracer_emissions(co_emis, 2)
    rate_parameters = build_musica_rate_parameters(co_emis, co_tracer_emis)
    times = np.arange(0.0, 50.0, 10.0)
    run_source_receptor(
        solver, state, rate_parameters, times, dt=10.0, K_diff=0.002, u_adv=0.006, dx=1.0
    )
