#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tagged tracer core functions for MUSICA/CheMPAS prototype.

Structure:
    tagged_tracer_core.py
        reusable parameterization routines (this module)

    examples/tagged_co_transport.py
        example driver script

    examples/tagged_co_transport.ipynb
        example notebook with plots

This first version handles:
    - one species: CO
    - one tagged tracer: CO_tracer
    - one tagged source grid cell
    - one receptor grid cell
    - total CO emissions everywhere
    - CO_tracer emissions only at the source cell
    - MUSICA local emissions/losses
    - explicit 1-D advection-diffusion transport after each MUSICA solve

Later, WRF/WRF-Chem/MPAS fields can replace the simple test inputs here.
"""

import numpy as np


def as_1d_array(value, n, name):
    """Return a 1-D array of length n from a scalar or length-n array."""
    arr = np.asarray(value, dtype=float)

    if arr.ndim == 0:
        return np.full(n, float(arr))

    if arr.shape != (n,):
        raise ValueError(f"{name} must be scalar or shape ({n},), got {arr.shape}")

    return arr


def build_1d_grid(num_grid_cells, source_cell, receptor_cell):
    """Build a simple 1-D grid and validate source/receptor indices."""
    grid = np.arange(num_grid_cells)

    if not 0 <= source_cell < num_grid_cells:
        raise ValueError("source_cell is outside the grid")

    if not 0 <= receptor_cell < num_grid_cells:
        raise ValueError("receptor_cell is outside the grid")

    return grid


def derive_simple_met(num_grid_cells, temperature=298.0, pressure=101325.0):
    """Create simple meteorological arrays for a MUSICA box-model test.

    This mirrors the dry-deposition example style, where met variables are
    defined/prepared before the parameterization is called.

    Returns
    -------
    dict
        temperature [K] and pressure [Pa].
    """
    return {
        "temperature": as_1d_array(temperature, num_grid_cells, "temperature"),
        "pressure": as_1d_array(pressure, num_grid_cells, "pressure"),
    }


def convert_surface_flux_to_concentration_tendency(surface_flux, layer_height):
    """Convert emission flux to concentration tendency.

    WRF/WRF-Chem-style placeholder:
        mol m-3 s-1 = (mol m-2 s-1) / layer_height

    Parameters
    ----------
    surface_flux : float or ndarray
        Emission flux in mol m-2 s-1.
    layer_height : float or ndarray
        Layer thickness in m.
    """
    surface_flux = np.asarray(surface_flux, dtype=float)
    layer_height = np.asarray(layer_height, dtype=float)

    if np.any(layer_height <= 0):
        raise ValueError("layer_height must be positive")

    return surface_flux / layer_height


def build_total_co_emissions(
    grid,
    base_emis=1.0e-12,
    enhancement=5.0e-13,
    enhancement_center=12.0,
    enhancement_width=30.0,
):
    """Build total CO emission tendency for every grid cell.

    Units: mol m-3 s-1.
    """
    grid = np.asarray(grid, dtype=float)

    return np.full(grid.shape, base_emis, dtype=float) + enhancement * np.exp(
        -((grid - enhancement_center) ** 2) / enhancement_width
    )


def build_single_cell_tracer_emissions(total_co_emis, source_cell):
    """Build CO_tracer emissions only at the tagged source cell."""
    total_co_emis = np.asarray(total_co_emis, dtype=float)

    tracer_emis = np.zeros_like(total_co_emis)
    tracer_emis[source_cell] = total_co_emis[source_cell]

    return tracer_emis


def build_musica_rate_parameters(co_emis, co_tracer_emis, co_loss=None, co_tracer_loss=None):
    """Build MUSICA user-defined rate parameter dictionary."""
    co_emis = np.asarray(co_emis, dtype=float)
    co_tracer_emis = np.asarray(co_tracer_emis, dtype=float)

    if co_loss is None:
        co_loss = np.zeros_like(co_emis)
    else:
        co_loss = np.asarray(co_loss, dtype=float)

    if co_tracer_loss is None:
        co_tracer_loss = np.zeros_like(co_tracer_emis)
    else:
        co_tracer_loss = np.asarray(co_tracer_loss, dtype=float)

    return {
        "EMIS.CO_emission": co_emis,
        "EMIS.CO_tracer_emission": co_tracer_emis,
        "LOSS.CO_loss": co_loss,
        "LOSS.CO_tracer_loss": co_tracer_loss,
    }


def transport_1d_advection_diffusion(C, K_diff, u_adv, dx, dt):
    """Apply simple 1-D advection-diffusion transport.

    This is a prototype transport operator, not a replacement for MPAS/WRF
    transport. It is included to prove the source-receptor tracer concept.
    """
    C = np.asarray(C, dtype=float)

    lap = np.zeros_like(C)
    lap[1:-1] = C[:-2] - 2.0 * C[1:-1] + C[2:]
    lap[0] = C[1] - C[0]
    lap[-1] = C[-2] - C[-1]

    diff_tendency = K_diff * lap / dx**2

    adv_tendency = np.zeros_like(C)

    if u_adv >= 0:
        adv_tendency[1:] = -u_adv * (C[1:] - C[:-1]) / dx
        adv_tendency[0] = 0.0
    else:
        adv_tendency[:-1] = -u_adv * (C[1:] - C[:-1]) / dx
        adv_tendency[-1] = 0.0

    C_new = C + dt * (diff_tendency + adv_tendency)

    return np.maximum(C_new, 0.0)


def source_receptor_proof_table(
    times, co_hist, co_tracer_hist, co_tracer_emis, source_cell, receptor_cell
):
    """Build a source-receptor proof table."""
    import pandas as pd

    return pd.DataFrame(
        {
            "time.s": times,
            "CO_at_receptor": co_hist[:, receptor_cell],
            "CO_tracer_at_source": co_tracer_hist[:, source_cell],
            "CO_tracer_at_receptor": co_tracer_hist[:, receptor_cell],
            "CO_tracer_emis_at_source": co_tracer_emis[source_cell],
            "CO_tracer_emis_at_receptor": co_tracer_emis[receptor_cell],
        }
    )


def history_to_dataframe(
    times,
    grid,
    co_hist,
    co_tracer_hist,
    co_emis,
    co_tracer_emis,
    temperature,
    pressure,
    source_cell,
    receptor_cell,
):
    """Convert history arrays to a tidy output DataFrame."""
    import pandas as pd

    rows = []

    for it, time_s in enumerate(times):
        for cell in grid:
            cell = int(cell)
            rows.append(
                {
                    "time.s": time_s,
                    "grid_cell": cell,
                    "is_source_cell": cell == source_cell,
                    "is_receptor_cell": cell == receptor_cell,
                    "CO.mol_m-3": co_hist[it, cell],
                    "CO_tracer.mol_m-3": co_tracer_hist[it, cell],
                    "CO_emis.mol_m-3_s-1": co_emis[cell],
                    "CO_tracer_emis.mol_m-3_s-1": co_tracer_emis[cell],
                    "temperature.K": temperature[cell],
                    "pressure.Pa": pressure[cell],
                }
            )

    return pd.DataFrame(rows)


def build_co_tracer_mechanism():
    """Build the MUSICA mechanism for CO and the tagged CO_tracer.

    The mechanism has an emission and a first-order loss for each species.
    The rate parameters ``EMIS.CO_emission``, ``EMIS.CO_tracer_emission``,
    ``LOSS.CO_loss``, and ``LOSS.CO_tracer_loss`` are user-defined.
    """
    import musica.mechanism_configuration as mc

    co = mc.Species(name="CO")
    co_tracer = mc.Species(name="CO_tracer")

    species = [co, co_tracer]
    gas = mc.Phase(name="gas", species=species)

    reactions = [
        mc.Emission(name="CO_emission", scaling_factor=1.0, products=[co], gas_phase=gas),
        mc.Emission(
            name="CO_tracer_emission", scaling_factor=1.0, products=[co_tracer], gas_phase=gas
        ),
        mc.FirstOrderLoss(name="CO_loss", scaling_factor=1.0, reactants=[co], gas_phase=gas),
        mc.FirstOrderLoss(
            name="CO_tracer_loss", scaling_factor=1.0, reactants=[co_tracer], gas_phase=gas
        ),
    ]

    return mc.Mechanism(
        name="single_co_tracer_transport",
        species=species,
        phases=[gas],
        reactions=reactions,
    )


def run_source_receptor(solver, state, rate_parameters, times, dt, K_diff, u_adv, dx):
    """Run the MUSICA solve and the explicit 1-D transport for each time step.

    The state must already have its conditions and initial concentrations.
    The function records the concentrations at the start of each step.

    Returns
    -------
    co_hist, co_tracer_hist : ndarray
        Concentrations in mol m-3, shape (len(times), num_grid_cells).
    """
    state.set_user_defined_rate_parameters(rate_parameters)

    co_hist = []
    co_tracer_hist = []

    for _ in times:
        conc_now = state.get_concentrations()
        co_hist.append(np.array(conc_now["CO"], dtype=float))
        co_tracer_hist.append(np.array(conc_now["CO_tracer"], dtype=float))

        # Local MUSICA box-model step
        solver.solve(state, dt)

        conc_after_musica = state.get_concentrations()

        # Explicit transport step
        state.set_concentrations(
            {
                name: transport_1d_advection_diffusion(
                    C=conc_after_musica[name], K_diff=K_diff, u_adv=u_adv, dx=dx, dt=dt
                )
                for name in ("CO", "CO_tracer")
            }
        )

    return np.asarray(co_hist), np.asarray(co_tracer_hist)
