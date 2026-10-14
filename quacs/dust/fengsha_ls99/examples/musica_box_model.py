"""MUSICA box model for the LS99 FENGSHA dust emission scheme.

The host driver gives an emission rate for each of 2 modes (accumulation, coarse) in g m-3 s-1.
The script multiplies by layer_height * jacobm (= layer_height / rjacm) and by 1e-3 to get a flux in kg m-2 s-1. Each flux is a MUSICA emission rate for one
species. MUSICA steps the species forward for 10 steps of 60 s.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.fengsha_ls99 import fengsha_dust_driver

N_BINS = 2
TITLE, MECHANISM_NAME, BIN_WORD = "LS99 FENGSHA", "quacs_dust_fengsha_ls99", "mode"
G_TO_KG = 1.0e-3

# One cell, source default inputs (dry, 80 % barren land, strong wind).
one = np.ones((1, 1))
DZ = 10.0  # layer thickness [m], also the driver's layer_height
JACOBM = 1.0  # host Jacobian field; the driver forms rjacm = 1 / jacobm
luf = np.zeros((1, 1, 24))  # full USGS24 land-use fractions; the driver picks the desert classes
luf[0, 0, 18] = 0.8  # USGS24 class 19: barren or sparsely vegetated, 80 % of the cell
out = fengsha_dust_driver(
    u10m=15.0 * one,
    soim_vol=0.02 * one,
    dens=1.2 * np.ones((1, 1, 1)),  # [x, level, y]
    snocov=0.0 * one,
    rn=0.0 * one,
    rc=0.0 * one,
    vegfrac=0.05 * one,
    soil_type=np.array([[1]]),
    luf=luf,
    jdate=2023100,  # YYYYDDD, day of year 100
    lwmask=one,
    layer_height=DZ * one,
    jacobm=JACOBM * one,
    scheme="usgs24",
)
# mode_mass is g m-3 s-1. Times layer_height * jacobm (= layer_height / rjacm) gives g m-2 s-1, times 1e-3 gives kg m-2 s-1.
F_BIN = G_TO_KG * out["mode_mass"][0, 0, :] * DZ * JACOBM  # kg m-2 s-1 per mode


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
