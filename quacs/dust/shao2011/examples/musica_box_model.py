"""MUSICA box model for the Shao (2011) dust emission scheme.

The scheme gives a dust flux for each of 5 size bins. Each flux is a MUSICA emission rate for
one species. MUSICA steps the species forward for 10 steps of 60 s.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.shao2011 import initialize, uoc_dust_driver
from quacs.dust.shao2011.shao2011 import N_SOIL

N_BINS = 5
TITLE, MECHANISM_NAME, BIN_WORD = "Shao 2011", "quacs_dust_shao2011", "bin"

# One cell, source default inputs (100 % sand, dry, strong wind), in the WRF layout.
one = np.ones((1, 1))
soil_top = np.zeros((1, N_SOIL, 1))  # [x, soil, y]: 100 % in soil class 1
soil_top[0, 0, 0] = 1.0
level = np.ones((1, 1, 1))
DT = 60.0
out = uoc_dust_driver(
    dust=np.zeros((1, 1, 1, N_BINS)),
    rho_phy=1.225 * level,
    dz8w=10.0 * level,
    smois=0.02 * level,
    ust=0.8 * one,
    isltyp=np.ones((1, 1), dtype=int),
    vegfra=5.0 * one,
    dt=DT,
    emis_dust=np.zeros((1, 1, N_BINS)),  # accumulated emitted mass [kg m-2], starts at zero
    soil_top_cat=soil_top,
    erod=np.array([[[0.5, 0.0, 0.0]]]),
    psd=initialize(),
)
# The driver adds flux * dt to emis_dust. Divide by dt to get the flux.
F_BIN = out["emis_dust"][0, 0, :] / DT  # kg m-2 s-1 per bin


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
