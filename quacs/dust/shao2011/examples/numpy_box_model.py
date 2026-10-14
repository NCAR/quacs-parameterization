"""NumPy box model for the Shao (2011) dust emission scheme.

One grid cell. The script prints the dry-threshold minimum, calls the host driver
`uoc_dust_driver`, and steps the dust mass in the lowest layer forward with dC = F / dz * dt.
The driver does not return the flux. It adds the emitted mass (flux * dt, kg m-2) to `emis_dust`.
The script starts `emis_dust` at zero and divides the result by dt to get the flux in kg m-2 s-1.
"""

import numpy as np

from quacs.dust.shao2011 import initialize, threshold_velocity_shao_lu, uoc_dust_driver
from quacs.dust.shao2011.shao2011 import DBIN, N_BINS, N_SOIL, RHOP_DEFAULT

# One cell, source default inputs (100 % sand, dry, strong wind).
ustar, gwet, cf, rhoa = 0.8, 0.02, 5.0, 1.225

# Dry-threshold curve: the particle diameter where the threshold is smallest.
d_curve = np.logspace(0, 3, 400)  # um
u_curve = threshold_velocity_shao_lu(d_curve, RHOP_DEFAULT / rhoa, rhoa)
d_min = d_curve[u_curve.argmin()]

# Host fields in WRF layout: [x, level, y] for the air, [x, soil, y] for the soil.
one = np.ones((1, 1))
level = np.ones((1, 1, 1))
dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
soil_top = np.zeros((1, N_SOIL, 1))  # soil class fractions, 100 % in class 1
soil_top[0, 0, 0] = 1.0
out = uoc_dust_driver(
    dust=np.zeros((1, 1, 1, N_BINS)),  # dust tracers [ug kg-1]
    rho_phy=rhoa * level,
    dz8w=dz * level,
    smois=gwet * level,
    ust=ustar * one,
    isltyp=np.ones((1, 1), dtype=int),
    vegfra=cf * one,
    dt=dt,
    emis_dust=np.zeros((1, 1, N_BINS)),  # accumulated emitted mass [kg m-2], starts at zero
    soil_top_cat=soil_top,
    erod=np.array([[[0.5, 0.0, 0.0]]]),
    psd=initialize(),
)
F_p = out["emis_dust"] / dt  # emitted mass / dt = kg m-2 s-1, shape (1, 1, 5)

print("Shao (2011) dust emission: box model")
print(f"  ustar = {ustar} m/s, gwet = {gwet}, cover = {cf} %, 100 % sand")
print(f"  smallest dry threshold: {u_curve.min():.3f} m/s at d = {d_min:.0f} um")
print(f"  threshold after cover and moisture: {out['ust_t_min'][0, 0]:.3f} m/s")
for n in range(N_BINS):
    lo = 0.0 if n == 0 else DBIN[n - 1]
    print(f"  dust flux, bin {n + 1} ({lo:.1f}-{DBIN[n]:.1f} um): {F_p[0, 0, n]:.4e} kg m-2 s-1")

mass = np.zeros(N_BINS)  # dust mass concentration per bin [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F_p[0, 0, :] / dz * dt
    print(f"{step:>5}  {F_p[0, 0, :].sum():>30.4e}  {mass.sum():>20.4e}")
