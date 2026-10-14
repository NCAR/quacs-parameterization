"""NumPy box model for the LS99 FENGSHA dust emission scheme.

One grid cell. The script calls the host driver `fengsha_dust_driver` for a dry sandy cell, shows
that snow and weak wind switch the emission off, and steps the dust mass in the lowest layer
forward with dC = F / dz * dt.

The driver returns `mode_mass`, the emission rate per mode in g m-3 s-1 (a rate per volume of the
lowest layer, weighted by the grid factor `rjacm` = 1 / `jacobm`). The script converts it to a flux in
kg m-2 s-1: F = mode_mass * layer_height * jacobm * 1e-3.
"""

import numpy as np

from quacs.dust.fengsha_ls99 import fengsha_dust_driver
from quacs.dust.fengsha_ls99.fengsha_ls99 import FRACMJ, FRACMK

G_TO_KG = 1.0e-3

dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]

# One cell, source default inputs (dry, 80 % barren land, strong wind).
one = np.ones((1, 1))
u10m = 15.0 * one
soim = 0.02 * one
dens = 1.2 * np.ones((1, 1, 1))  # [x, level, y]
snow = 0.0 * one
rain = 0.0 * one
veg = 0.05 * one
stype = np.array([[1]])
luf = np.zeros((1, 1, 24))  # full USGS24 land-use fractions; the driver picks the desert classes
luf[0, 0, 18] = 0.8  # USGS24 class 19: barren or sparsely vegetated, 80 % of the cell


def flux(**kw):
    args = dict(
        u10m=u10m,
        soim_vol=soim,
        dens=dens,
        snocov=snow,
        rn=rain,
        rc=0.0 * one,
        vegfrac=veg,
        soil_type=stype,
        luf=luf,
        jdate=2023100,  # YYYYDDD, day of year 100
        lwmask=one,  # land
        layer_height=dz * one,
        jacobm=one,
        scheme="usgs24",
    )
    args.update(kw)
    out = fengsha_dust_driver(**args)
    # mode_mass: g m-3 s-1. Times layer_height * jacobm (= layer_height / rjacm) gives g m-2 s-1.
    return out["mode_mass"] * args["layer_height"][..., None] * args["jacobm"][..., None]


F_g = flux()  # g m-2 s-1, shape (1, 1, 2)
print("LS99 FENGSHA dust emission: box model")
print(f"  u10 = {u10m[0, 0]} m/s, soil moisture = {soim[0, 0]}, barren = 80 %")
print(f"  mode split, accumulation / coarse: {FRACMJ} / {FRACMK}")
print(f"  total dust flux, dry cell:      {F_g.sum():.4e} g m-2 s-1")
print(f"  total dust flux, snow 0.5:      {flux(snocov=0.5 * one).sum():.4e} g m-2 s-1")
print(f"  total dust flux, wind 3 m/s:    {flux(u10m=3.0 * one).sum():.4e} g m-2 s-1")

F_p = F_g * G_TO_KG  # kg m-2 s-1
mass = np.zeros(2)  # dust mass concentration per mode [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F_p[0, 0, :] / dz * dt
    print(f"{step:>5}  {F_p[0, 0, :].sum():>30.4e}  {mass.sum():>20.4e}")
