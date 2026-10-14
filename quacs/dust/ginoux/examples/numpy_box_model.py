"""NumPy box model for the GOCART-Ginoux dust emission scheme.

One grid cell. The script prints the bin thresholds, calls the host driver
`ginoux_dust_driver`, and steps the dust mass in the lowest layer forward with dC = F / dz * dt.
The driver gets the soil moisture as a volume fraction and divides by the soil porosity to get
`gwet`. It returns the flux in kg m-2 s-1, so no unit change is needed.
"""

import numpy as np

from quacs.dust.ginoux import ginoux_dust_driver, threshold_velocity
from quacs.dust.ginoux.ginoux import DEN_DUST, FRAC_S, LAND, NDCLS, POROSITY, REFF_DUST

# One cell, source default inputs.
w10m, gwet, erod, rhoa = 12.0, 0.1, 0.5, 1.2

u_t = threshold_velocity(REFF_DUST, DEN_DUST, np.array([[rhoa]]))

print("Ginoux (GOCART) dust emission: box model")
print(f"  w10m = {w10m} m/s, gwet = {gwet}, erodibility = {erod}")
print(f"  bin radius [um]: {REFF_DUST * 1e6}")
print(f"  bin source fraction [-]: {FRAC_S} (sum = {FRAC_S.sum():.2f})")
print(f"  threshold wind speed per bin [m/s]: {u_t[0, 0, :].round(3)}")

# Host fields in WRF layout: [x, level, y] for the air, [x, soil layer, y] for the soil.
one = np.ones((1, 1))
level = np.ones((1, 1, 1))
dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
isltyp = 1  # USGS soil type 1 (sand)
out = ginoux_dust_driver(
    dust=np.zeros((1, 1, 1, len(FRAC_S))),  # dust tracers [ug kg-1]
    rho_phy=rhoa * level,
    dz8w=dz * level,  # below 12 m, so the driver uses the model wind u_phy, v_phy
    dt=dt,
    u10=w10m * one,
    v10=0.0 * one,
    u_phy=w10m * level,
    v_phy=0.0 * level,
    smois=gwet * POROSITY[isltyp - 1] * level,  # volume fraction; gwet = smois / porosity
    isltyp=isltyp * one.astype(int),
    oro=np.full((1, 1), LAND),
    erod=np.full((1, 1, NDCLS), erod),
    emis_dust=np.zeros((1, 1, len(FRAC_S))),
    dustin=np.zeros((1, 1, len(FRAC_S))),
)
F_p = out["flux"]  # kg m-2 s-1, shape (1, 1, 5)
print(f"  dust flux per bin [kg m-2 s-1]: {F_p[0, 0, :]}")

mass = np.zeros(len(FRAC_S))  # dust mass concentration per bin [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F_p[0, 0, :] / dz * dt
    print(f"{step:>5}  {F_p[0, 0, :].sum():>30.4e}  {mass.sum():>20.4e}")
