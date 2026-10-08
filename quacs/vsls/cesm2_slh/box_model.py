"""
MUSICA box model for the CESM2-SLH sea-salt bromine reactions.

Python computes the rate constants with :func:`seasalt_rate_constants`.
MICM integrates the concentrations.
"""

from typing import Dict

import musica
import numpy as np

from .mechanism import build_mechanism, rate_parameter_key
from .seasalt import R_GAS, seasalt_rate_constants


def run_box_model(
    met: Dict[str, float],
    initial_pptv: Dict[str, float],
    time_step_s: float = 600.0,
    n_steps: int = 144,
) -> Dict[str, object]:
    """Run the sea-salt bromine mechanism in one cell.

    Parameters
    ----------
    met : dict
        Keyword arguments for :func:`seasalt_rate_constants` (scalars).
    initial_pptv : dict
        Initial mixing ratio (pptv) for each species. Missing species start at 0.
    time_step_s : float
        MICM time step (s).
    n_steps : int
        Number of time steps.

    Returns
    -------
    dict
        ``rates`` (reaction name -> k, 1/s), ``n_air`` (mol/m3),
        ``times`` (s, shape ``(n_steps + 1,)``), and ``concentrations``
        (species name -> mol/m3, shape ``(n_steps + 1,)``).
    """
    rates = {name: float(k) for name, k in seasalt_rate_constants(**met).items()}

    mechanism, species = build_mechanism()
    solver = musica.MICM(
        mechanism=mechanism, solver_type=musica.SolverType.rosenbrock_standard_order
    )
    state = solver.create_state(1)
    state.set_conditions(temperatures=[met["temperature_k"]], pressures=[met["pressure_pa"]])

    n_air = met["pressure_pa"] / (R_GAS * met["temperature_k"])  # mol/m3
    c0 = {s: initial_pptv.get(s, 0.0) * 1e-12 * n_air for s in species}
    state.set_concentrations({s: [v] for s, v in c0.items()})
    state.set_user_defined_rate_parameters(
        {rate_parameter_key(name): [k] for name, k in rates.items()}
    )

    history = {s: [v] for s, v in c0.items()}
    for _ in range(n_steps):
        solver.solve(state, time_step_s)
        conc = state.get_concentrations()
        for s in history:
            history[s].append(conc[s][0])

    return {
        "rates": rates,
        "n_air": n_air,
        "times": np.arange(n_steps + 1) * time_step_s,
        "concentrations": {s: np.asarray(v) for s, v in history.items()},
    }
