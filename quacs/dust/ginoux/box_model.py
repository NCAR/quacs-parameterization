"""
MUSICA box model for the GOCART dust emission scheme (Ginoux et al., 2001).

The physics is in ginoux.py. This module builds the MICM mechanism, makes
the solver, and steps it in time. It gets F_p from the scheme, divides it
by the layer thickness dz to get a volumetric emission rate, and sets that
rate on the MICM state with the EMIS.* user-defined parameters:

    EMIS.dust_emis_bin_i = F_p[i] / dz      [kg m-3 s-1]

The dust bins are tracers with a placeholder molecular weight of 1 kg/mol,
so the MICM concentration has the same value as the dust mass in kg m-3.
"""

import numpy as np

from .ginoux import (
    C_DEFAULT,
    LAND,
    RADIUS_DEFAULT,
    RHOP_DEFAULT,
    SP_DEFAULT,
    dust_emission,
    threshold_velocity,
)

DUST_BIN_MW = 1.0  # kg/mol, placeholder so that mol m-3 == kg m-3


def build_mechanism(n_bins=len(SP_DEFAULT)):
    """Return a MICM mechanism with one dust tracer and one emission for each bin."""
    import musica.mechanism_configuration as mc

    dust_species = [
        mc.Species(name=f"dust_bin_{i + 1}", molecular_weight_kg_mol=DUST_BIN_MW)
        for i in range(n_bins)
    ]
    gas = mc.Phase(name="gas", species=dust_species)
    emissions = [
        mc.Emission(
            name=f"dust_emis_bin_{i + 1}",
            scaling_factor=1.0,
            products=[dust_species[i]],
            gas_phase=gas,
        )
        for i in range(n_bins)
    ]
    return mc.Mechanism(
        name="quacs_dust_ginoux",
        species=dust_species,
        phases=[gas],
        reactions=emissions,
    )


def emission_rates(w10m, u_t, gwet, oro, frlake, S, s_p, dz, C=C_DEFAULT):
    """Convert the surface flux F_p [kg m-2 s-1] to MICM emission rates [kg m-3 s-1]."""
    F_p = dust_emission(w10m=w10m, u_t=u_t, gwet=gwet, oro=oro, frlake=frlake, S=S, s_p=s_p, C=C)
    return {f"EMIS.dust_emis_bin_{i + 1}": [F_p[0, 0, i] / dz] for i in range(len(s_p))}


def _cell(value):
    return np.array([[float(value)]])


def run_box_model(
    w10m=12.0,
    gwet=0.1,
    S=0.5,
    frlake=0.0,
    rhoa=1.2,
    dz=50.0,
    dt=120.0,
    n_steps=30,
    temperature=300.0,
    pressure=101325.0,
    C=C_DEFAULT,
):
    """Run the dust emission box model for one land cell.

    Returns
    -------
    dict with keys
        ``times`` : (n_steps + 1,) time [s]
        ``concs`` : (n_steps + 1, n_bins) dust mass for each bin [kg m-3]
        ``rates`` : (n_bins,) emission rate for each bin [kg m-3 s-1]
        ``u_t``   : bin-mean dry-soil threshold wind speed [m s-1]
    """
    import musica

    s_p = SP_DEFAULT
    n_bins = len(s_p)

    # Dry-soil threshold for each bin (Marticorena & Bergametti, 1995).
    # dust_emission takes one u_t for each cell, so use the bin mean.
    u_t_bins = threshold_velocity(RADIUS_DEFAULT, RHOP_DEFAULT, _cell(rhoa))
    u_t = _cell(u_t_bins[0, 0, :].mean())

    rates = emission_rates(
        w10m=_cell(w10m),
        u_t=u_t,
        gwet=_cell(gwet),
        oro=_cell(LAND),
        frlake=_cell(frlake),
        S=_cell(S),
        s_p=s_p,
        dz=dz,
        C=C,
    )

    solver = musica.MICM(
        mechanism=build_mechanism(n_bins),
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(number_of_grid_cells=1)
    state.set_conditions([temperature], [pressure])
    state.set_concentrations({f"dust_bin_{i + 1}": [0.0] for i in range(n_bins)})
    state.set_user_defined_rate_parameters(rates)

    times = np.arange(n_steps + 1) * dt
    concs = np.zeros((n_steps + 1, n_bins))
    for step in range(n_steps):
        solver.solve(state, dt)
        c = state.get_concentrations()
        concs[step + 1] = [c[f"dust_bin_{i + 1}"][0] for i in range(n_bins)]

    return {
        "times": times,
        "concs": concs,
        "rates": np.array([rates[f"EMIS.dust_emis_bin_{i + 1}"][0] for i in range(n_bins)]),
        "u_t": float(u_t[0, 0]),
    }
