"""
Run the GOCART dust emission scheme in a MUSICA box model for one land cell.

The emission rate is constant, so the dust mass must grow linearly in time.
The script prints the MICM result next to that closed-form value.
"""

from quacs.dust.ginoux import C_DEFAULT, run_box_model

w10m = 12.0  # 10-m wind speed [m s-1]
gwet = 0.1  # gravimetric soil moisture [-]
dz = 50.0  # lowest-layer thickness [m]
dt = 120.0  # time step [s]
n_steps = 30

out = run_box_model(w10m=w10m, gwet=gwet, dz=dz, dt=dt, n_steps=n_steps)

print("--- MUSICA dust-emission box model (GOCART) ---")
print(f"  w10m={w10m:.1f} m/s, gwet={gwet:.2f}, u_t (bin-mean)={out['u_t']:.3f} m/s")
print(f"  dz={dz:.1f} m, dt={dt:.0f} s, n_steps={n_steps}, C={C_DEFAULT}\n")
print(f"{'time (min)':>10}  {'total dust [kg/m3]':>22}  {'expected [kg/m3]':>22}")
print("-" * 60)

expected = out["rates"].sum() * out["times"]
for t, actual, exp in zip(out["times"], out["concs"].sum(axis=1), expected):
    print(f"{t / 60:>10.2f}  {actual:>22.6e}  {exp:>22.6e}")
