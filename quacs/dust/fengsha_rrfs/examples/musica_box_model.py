"""MUSICA box model for the FENGSHA RRFS dust emission scheme.

The scheme gives a dust flux for each of 5 size bins. Each flux is a MUSICA emission rate for
one species. MUSICA steps the species forward for 10 steps of 60 s.
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.fengsha_rrfs import rrfs_dust_step
from quacs.dust.fengsha_rrfs.fengsha_rrfs import CON_RD

N_BINS = 5
TITLE, MECHANISM_NAME, BIN_WORD = "FENGSHA RRFS", "quacs_dust_fengsha_rrfs", "bin"

# One cell with example inputs. The source has no default for these; the host model supplies them.
one = np.ones((1, 1))
T_AIR = 300.0
F_BIN = rrfs_dust_step(
    dust12m=np.tile([0.15, 0.60, 0.55, 0.80, 0.30], (1, 1, 12, 1)),  # clay, rdrag, sand, ssm, uthrs
    month=1,
    zorl=np.array([[1.0]]),  # cm
    snow=np.array([[0.0]]),  # mm
    land=np.array([[1]]),
    smc=np.array([[[0.3, 0.02]]]),  # layer 2 is used
    tslb=np.array([[[290.0]]]),  # layer 1 is used
    pr3d=np.array([[[1.0e5, 0.99e5]]]),
    prl3d=np.array([[[1.2 * CON_RD * T_AIR]]]),  # rho_a = p / (Rd T) = 1.2 kg m-3
    tk3d=np.array([[[T_AIR]]]),
    ustar=1.0 * one,
    garea=one,
    soiltyp=np.array([[2]]),
    dt=1.0,
    dust=np.zeros((1, 1)),  # RRFS dust tracer, ug/kg
    coarse_pm=np.zeros((1, 1)),  # RRFS coarse PM tracer, ug/kg
).rate[0, 0, :]  # kg m-2 s-1 per bin


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
