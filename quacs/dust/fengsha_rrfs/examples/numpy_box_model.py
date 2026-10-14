"""NumPy box model for the FENGSHA RRFS dust emission scheme.

One grid cell. The script prints the threshold, moisture factor and bin fractions, calls the
scheme, and steps the dust mass in the lowest layer forward with dC = F / dz * dt.
"""

import numpy as np

from quacs.dust.fengsha_rrfs import (
    RrfsSettings,
    flux_ratio_mb95,
    kok_bin_fractions,
    moisture_correction_fecan,
    rrfs_dust_step,
)
from quacs.dust.fengsha_rrfs.fengsha_rrfs import CON_RD, N_BINS, REFF_DUST

# One cell with example inputs. The source has no default for these; the host model supplies
# them. With ustar 0.9 m/s the drag-reduced wind (0.54 m/s) stays below the wet threshold
# (0.5475 m/s) and the flux is zero, so the example uses 1.0 m/s.
ustar, slc, rho_a, clay, sand = 1.0, 0.02, 1.2, 0.15, 0.55
ssm, rdrag, uthrs = 0.80, 0.60, 0.30

settings = RrfsSettings()
slc_eff = settings.moist_corr * slc  # the scheme scales the soil moisture before the correction
h = moisture_correction_fecan(np.array([[slc_eff]]), np.array([[sand]]), np.array([[clay]]))[0, 0]
kvh = flux_ratio_mb95(np.array([[clay]]))[0, 0]
kok = kok_bin_fractions()

# The driver gets p and T and derives the air density. p is set so that rho_a = p / (Rd T) is 1.2.
T_air = 300.0
one = np.ones((1, 1))
# rrfs_dust_step takes host arrays: 12 months of soil fields, layers on the last axis, and the two
# RRFS dust tracers (dust and coarse PM, ug/kg). The wrapper builds the 5 bins from them.
dust12m = np.tile([clay, rdrag, sand, ssm, uthrs], (1, 1, 12, 1))  # W:818-822 field order
result = rrfs_dust_step(
    dust12m=dust12m,
    month=1,
    zorl=np.array([[1.0]]),  # cm, the step converts to 0.01 m (W:814)
    snow=np.array([[0.0]]),  # mm
    land=np.array([[1]]),
    smc=np.array([[[0.3, slc]]]),  # soil layers; layer 2 is used
    tslb=np.array([[[290.0]]]),  # soil layers; layer 1 is used
    pr3d=np.array([[[1.0e5, 0.99e5]]]),  # layer edges
    prl3d=np.array([[[rho_a * CON_RD * T_air]]]),  # layer middles; sets rho_a = p / (Rd T)
    tk3d=np.array([[[T_air]]]),
    ustar=np.array([[ustar]]),
    garea=one,
    soiltyp=np.array([[2]]),
    dt=1.0,
    dust=np.zeros((1, 1)),  # RRFS dust tracer, ug/kg
    coarse_pm=np.zeros((1, 1)),  # RRFS coarse PM tracer, ug/kg
)
F_p = result.rate  # kg m-2 s-1, shape (1, 1, 5)

print("FENGSHA RRFS dust emission: box model")
print(f"  ustar = {ustar} m/s, soil moisture = {slc}, clay = {clay}, sand = {sand}")
print(f"  dry threshold = {uthrs} m/s, drag partition = {rdrag}")
print(f"  moisture factor (Fecan) = {h:.4f}, wet threshold = {uthrs * h:.4f} m/s")
print(f"  vertical-to-horizontal flux ratio = {kvh:.4e} m-1")
print(f"  bin radius [um]: {REFF_DUST * 1e6}")
print(f"  bin fraction [-]: {np.round(kok, 4)} (sum = {np.sum(kok):.4f})")
print(f"  dust flux per bin [kg m-2 s-1]: {F_p[0, 0, :]}")

dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
mass = np.zeros(N_BINS)  # dust mass concentration per bin [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F_p[0, 0, :] / dz * dt
    print(f"{step:>5}  {F_p[0, 0, :].sum():>30.4e}  {mass.sum():>20.4e}")
