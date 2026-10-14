"""MUSICA box model for the Leung (2023) dust emission scheme.

The driver gives a dust flux for each of 4 size bins. Each flux is a MUSICA emission rate for
one species. MUSICA steps the species forward for 10 steps of 60 s.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.leung2023 import leung_dust_driver

N_BINS = 4
TITLE, MECHANISM_NAME, BIN_WORD = "Leung 2023", "quacs_dust_leung2023", "bin"

# One cell, inputs of the README usage example.
one = np.ones((1, 1))
result = leung_dust_driver(
    ustar=0.7 * one,
    rhoa=1.1 * one,
    h2osoi_vol=0.02 * one,
    tlai=0.06 * one,
    tsai=0.04 * one,
    frac_sno=0.0 * one,
    pct_sand=65.0 * one,
    pct_clay=10.0 * one,
    organic=0.0 * one,
    h2osoi_liq=1.0 * one,
    h2osoi_ice=0.0 * one,
    land_fraction=one,
    landunit_type=np.array([[1]]),
    patch_active=np.array([[True]]),
    patch_wtlunit=one,
    rockfrc=0.5 * one,
    vegefrc=0.1 * one,
    sparfrc=0.1 * one,
    roughfct=0.7 * one,
    forc_pbot=1.0e5 * one,
    forc_t=295.0 * one,
    ram1=100.0 * one,
)
F_BIN = result["flx_mss_vrt_dst"][0, 0, :]  # kg m-2 s-1 per bin


# MUSICA works in mol m-3. The molecular weight below is a nominal value for this box only.
# It is not a property of dust. The mass results do not depend on it.
MW = 0.1  # kg mol-1, nominal

species = [mc.Species(name=f"dust_bin_{i + 1}", molecular_weight_kg_mol=MW) for i in range(N_BINS)]
gas = mc.Phase(name="gas", species=species)
emissions = [
    mc.Emission(
        name=f"dust_emis_bin_{i + 1}", scaling_factor=1.0, products=[species[i]], gas_phase=gas
    )
    for i in range(N_BINS)
]
mechanism = mc.Mechanism(name=MECHANISM_NAME, species=species, phases=[gas], reactions=emissions)

solver = musica.MICM(mechanism=mechanism, solver_type=musica.SolverType.rosenbrock_standard_order)
state = solver.create_state(number_of_grid_cells=1)
state.set_conditions([300.0], [101325.0])  # temperature [K], pressure [Pa]
state.set_concentrations({sp.name: [0.0] for sp in species})

dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
rates = {f"EMIS.dust_emis_bin_{i + 1}": [F_BIN[i] / dz / MW] for i in range(N_BINS)}  # mol m-3 s-1
state.set_user_defined_rate_parameters(rates)

print(f"{TITLE}: MUSICA box model, {N_BINS} dust {BIN_WORD}s")
print(f"  total dust flux = {F_BIN.sum():.4e} kg m-2 s-1, nominal MW = {MW} kg/mol")
for step in range(1, 11):
    solver.solve(state, dt)
    conc = state.get_concentrations()
    mass = [conc[f"dust_bin_{i + 1}"][0] * MW for i in range(N_BINS)]  # kg m-3
    cols = "  ".join(f"{m:.3e}" for m in mass)
    print(f"step {step:2d}  t = {step * dt:4.0f} s  dust per {BIN_WORD} [kg m-3]: {cols}")
