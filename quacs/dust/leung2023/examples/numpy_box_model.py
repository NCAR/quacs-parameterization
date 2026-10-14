"""NumPy box model for the Leung (2023) dust emission scheme.

One grid cell. The script calls `leung_dust_driver`, prints the threshold, drag, moisture and
intermittency values from its diagnostics, and steps the dust mass in the lowest layer forward
with dC = F / dz * dt. F is the total of the 4 size bins, the same flux as the MUSICA example.
The 4 bins hold particles up to 10 um, so F is less than `flux_before_split` (the Fortran source
drops the larger particles).
"""

import numpy as np

from quacs.dust.leung2023 import leung_dust_driver

# One cell, inputs of the README usage example.
one = np.ones((1, 1))
ustar, h2osoi_vol, tlai, tsai, clay = 0.7, 0.02, 0.06, 0.04, 10.0
result = leung_dust_driver(
    ustar=ustar * one,
    rhoa=1.1 * one,
    h2osoi_vol=h2osoi_vol * one,
    tlai=tlai * one,
    tsai=tsai * one,
    frac_sno=0.0 * one,
    pct_sand=65.0 * one,
    pct_clay=clay * one,
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
F = float(result["flx_mss_vrt_dst_tot"][0, 0])  # kg m-2 s-1, sum of the 4 bins
d = result

print("Leung (2023) dust emission: box model")
print(
    f"  ustar = {ustar} m/s, soil moisture = {h2osoi_vol}, LAI+SAI = {tlai + tsai}, clay = {clay} %"
)
print(f"  dry fluid threshold u*ft0 = {d['u_ft0'][0, 0]:.4f} m/s")
print(
    f"  wet fluid threshold u*ft  = {d['u_ft'][0, 0]:.4f} m/s (moisture factor {d['fm'][0, 0]:.3f})"
)
print(f"  impact threshold u*it     = {d['u_it'][0, 0]:.4f} m/s")
print(f"  ustar after drag u*s      = {d['u_s'][0, 0]:.4f} m/s")
print(f"  standardized u*st         = {d['u_st'][0, 0]:.4f} m/s")
print(f"  emission coefficient Cd   = {d['Cd'][0, 0]:.3e}")
print(f"  intermittency eta         = {d['eta'][0, 0]:.3f}")
print(f"  flux before the bin split = {float(d['flux_before_split'][0, 0]):.4e} kg m-2 s-1")
print(f"  total dust flux, 4 bins   = {F:.4e} kg m-2 s-1")

dz, dt = 10.0, 60.0  # layer thickness [m], time step [s]
mass = 0.0  # dust mass concentration [kg m-3]
print(f"{'step':>5}  {'total dust flux [kg m-2 s-1]':>30}  {'dust mass [kg m-3]':>20}")
for step in range(1, 11):
    mass += F / dz * dt
    print(f"{step:>5}  {F:>30.4e}  {mass:>20.4e}")
