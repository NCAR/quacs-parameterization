"""NumPy box model for the AFWA dust emission scheme.

One grid cell. The script prints the saltation thresholds, calls the host driver
`afwa_dust_driver`, and steps the dust mass in the lowest layer forward with dC = F / dz * dt.
The driver returns the flux in kg m-2 s-1, so no unit change is needed.
"""

import numpy as np

from quacs.dust.afwa import (
    afwa_dust_driver,
    saltation_area_fraction,
    threshold_friction_velocity,
)
from quacs.dust.afwa.afwa import DEN_SALT, DISTR_DUST, LAND, NDUST, REFF_SALT

# One cell, source default inputs (dry sandy soil, strong wind).
ustar, volsm, airden, erod, clay, sand = 0.6, 0.02, 1.225, 0.8, 0.20, 0.60

u_ts0 = threshold_friction_velocity(REFF_SALT, DEN_SALT, np.array([[airden]]))[0, 0, :]
ds_rel = saltation_area_fraction(np.array([[clay]]), np.array([[sand]]))[0, 0, :]

print("AFWA dust emission: box model")
print(f"  ustar = {ustar} m/s, soil moisture = {volsm}, erodibility = {erod}")
print(f"  saltation bin radius [um]: {REFF_SALT * 1e6}")
print(f"  dry threshold friction velocity [m/s]: {u_ts0.round(3)}")
print(f"  saltation area fraction [-]: {ds_rel.round(4)} (sum = {ds_rel.sum():.3f})")
print(f"  dust bin split [-]: {DISTR_DUST} (sum = {DISTR_DUST.sum():.4f})")

# Host fields in WRF layout: [x, level, y] for the air, [x, soil layer, y] for the soil.
one = np.ones((1, 1))
level = np.ones((1, 1, 1))
dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
out = afwa_dust_driver(
    dust=np.zeros((1, 1, 1, NDUST)),  # dust tracers [ug kg-1]
    rho_phy=airden * level,
    dz8w=dz * level,
    smois=volsm * level,
    dt=dt,
    u10=8.0 * one,  # 10-m wind [m s-1], used only for the lofting diagnostic
    v10=0.0 * one,
    emis_dust=np.zeros((1, 1, NDUST)),
    tot_edust=np.zeros((1, 1)),
    dustin=np.zeros((1, 1, NDUST)),
    ustar=ustar * one,
    znt=0.01 * one,
    snowh=0.0 * one,
    erod=np.array([[[erod, 0.0, 0.0]]]),
    clay=clay * one,
    sand=sand * one,
    oro=np.full((1, 1), LAND),
    isltyp=3 * one,
)
F_p = out["flux"]  # kg m-2 s-1, shape (1, 1, 5)
print(f"  lofting potential = {out['afwa_dustloft'][0, 0]:.3f} m/s")
print(f"  dust flux per bin [kg m-2 s-1]: {F_p[0, 0, :]}")

mass = np.zeros(NDUST)  # dust mass concentration per bin [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F_p[0, 0, :] / dz * dt
    print(f"{step:>5}  {F_p[0, 0, :].sum():>30.4e}  {mass.sum():>20.4e}")
