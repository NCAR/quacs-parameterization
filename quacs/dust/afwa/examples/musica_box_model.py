"""MUSICA box model for the AFWA dust emission scheme.

The scheme gives a dust flux for each of 5 size bins. Each flux is a MUSICA emission rate for
one species. MUSICA steps the species forward for 10 steps of 60 s.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.afwa import afwa_dust_driver

N_BINS = 5
TITLE, MECHANISM_NAME, BIN_WORD = "AFWA", "quacs_dust_afwa", "bin"

# One cell, source default inputs, in the WRF layout the host driver takes.
one = np.ones((1, 1))
level = np.ones((1, 1, 1))
out = afwa_dust_driver(
    dust=np.zeros((1, 1, 1, N_BINS)),
    rho_phy=1.225 * level,
    dz8w=10.0 * level,
    smois=0.02 * level,
    dt=60.0,
    u10=8.0 * one,
    v10=0.0 * one,
    emis_dust=np.zeros((1, 1, N_BINS)),
    tot_edust=np.zeros((1, 1)),
    dustin=np.zeros((1, 1, N_BINS)),
    ustar=0.6 * one,
    znt=0.01 * one,
    snowh=0.0 * one,
    erod=np.array([[[0.8, 0.0, 0.0]]]),
    clay=0.20 * one,
    sand=0.60 * one,
    oro=one,
    isltyp=3 * one,
)
F_BIN = out["flux"][0, 0, :]  # kg m-2 s-1 per bin


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
