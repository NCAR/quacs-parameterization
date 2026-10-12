# -*- coding: utf-8 -*-
"""
musica_box_model.py
===================
MUSICA box-model for the Leung et al. (2023/2024) dust emission scheme
(CESM2/CLM5 implementation).

The physics lives in quacs/dust/leung/leung.py.

This file builds the MICM mechanism, creates the solver, and time-steps.
Each step it asks the science module for F_p [x, y, n_bins],
divides by dz (and the placeholder molecular weight) to get volumetric
molar emission rates, and sets them on the MICM state via the EMIS.*
user-defined parameter convention.

LOSS:  state.set_user_defined_rate_parameters({"LOSS.Hg0_drydep": [k]})
EMIS:  state.set_user_defined_rate_parameters({"EMIS.dust_emis_bin_i": [e_i]})

Units
-----
    F_p                [kg m-2 s-1]   surface flux per bin (science module)
    F_p / dz           [kg m-3 s-1]   volumetric mass source
    F_p / dz / MW      [mol m-3 s-1]  what MICM's EMIS rate expects
    MICM concentration [mol m-3]      -> x MW -> [kg m-3] for output

Key differences from Ginoux / Zender / K14 box models
------------------------------------------------------
The Leung et al. scheme is the most physically complete of the four:

  1. 4 transport bins (0.1-1, 1-2.5, 2.5-5, 5-10 µm), same as Zender.
  2. Shao & Lu (2000) dry threshold at Dp = 127 µm (vs 75 µm in Zender/K14).
  3. Explicit drag partition: rock (Marticorena & Bergametti 1995) +
     vegetation (Okin 2008) reduce u* at the soil surface via satellite-
     derived z0a and VAI-based scheme.
  4. Impact (dynamic) threshold u*it = Bit * u*ft0 (Comola et al. 2019)
     replaces the fluid threshold in the emission equation, enabling
     intermittent emission in marginal regions.
  5. Intermittency factor eta ∈ [0,1] from a Gaussian turbulence model
     (Comola et al. 2019), driven by the Obukhov length obu.
  6. Upscaling correction map Kc_map (Leung et al. 2024 Eq. 12) applied
     as a time-invariant spatial scaling.
  7. Extended VAI threshold = 1.0 (vs 0.3 in Z03/K14).  Soil-moisture
     threshold and clay term follow CTSM: gwc_thr = 0.17 + 0.14*fclay
     (dust_moist_fact = 1) and fclay' = 0.1 + 0.5*min(fclay, 0.2).

The static and synthetic-dynamic inputs below mirror the standalone demo in
examples/numpy_box_model.py: same numbers, same diagnostics, solved by the Rosenbrock
integrator instead of an explicit Forward Euler step.

References:
    Leung, D. M., et al. (2023). A new process-based and scale-aware desert
    dust emission scheme – Part I. Atmos. Chem. Phys., 23, 6487–6523.
    https://doi.org/10.5194/acp-23-6487-2023

    Leung, D. M., et al. (2024). A new process-based and scale-aware desert
    dust emission scheme – Part II. Atmos. Chem. Phys., 24, 2287–2318.
    https://doi.org/10.5194/acp-24-2287-2024

    Comola, F., et al. (2019). The Intermittency of Wind-Driven Sand Transport.
    Geophys. Res. Lett., 46, 13430–13440.
    https://doi.org/10.1029/2019GL085739
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.leung.leung import (
    BIT,
    CD0,
    CE,
    CK,
    DP_MED,
    DST_SRC_NBR,
    KAPPA_MAX,
    MSS_FRC_SRC,
    NDST,
    OVR_SRC_SNK_MSS,
    RHO_AIR0,
    U_FT0_DRY_REF,
    U_ST0,
    _log_law_factor,
    clay_mass_fraction_leung,
    drag_partition,
    dust_emission,
    intermittency_factor,
    moisture_threshold_factor,
    sl2000_u_ft0_dry,
    turbulent_sigma_us,
)

# ---------------------------------------------------------------------------
# 1. Scheme metadata
# ---------------------------------------------------------------------------
n_bins = NDST  # 4 transport bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm
print(f"Leung et al. 2023/2024 scheme: {n_bins} transport bins")

# ---------------------------------------------------------------------------
# 2. MICM mechanism: one emitted species per dust bin
# ---------------------------------------------------------------------------
# Molecular weight is a dimensional placeholder (kg/mol).
# MICM works in mol m-3 (concentrations) and mol m-3 s-1 (EMIS rates).
# All kg <-> mol conversions below divide/multiply by DUST_BIN_MW explicitly,
# so results stay correct if this placeholder is ever changed.
DUST_BIN_MW = 1.0  # kg/mol, dimensional placeholder

dust_species = [
    mc.Species(name=f"dust_bin_{i + 1}", molecular_weight_kg_mol=DUST_BIN_MW) for i in range(n_bins)
]
gas = mc.Phase(name="gas", species=dust_species)

dust_emissions = [
    mc.Emission(
        name=f"dust_emis_bin_{i + 1}",
        scaling_factor=1.0,
        products=[dust_species[i]],
        gas_phase=gas,
    )
    for i in range(n_bins)
]

mechanism = mc.Mechanism(
    name="quacs_dust_leung2024",
    species=dust_species,
    phases=[gas],
    reactions=dust_emissions,
)

# ---------------------------------------------------------------------------
# 3. Static inputs (single cell [1, 1])
#    Mirrors the demo in examples/numpy_box_model.py.
# ---------------------------------------------------------------------------
# Soil properties (static)
fclay = np.array([[0.10]])  # clay mass fraction [-], 10 %
z0a = np.array([[1.0e-4]])  # aeolian roughness length [m] (Prigent 2005)
frac_rock = np.array([[0.0]])  # fractional rock area Ar [-], no rocks
frac_veg = np.array([[0.0]])  # fractional vegetation area Av [-], no plants
Kc_map = np.array([[1.0]])  # upscaling correction factor (unity = no correction)
is_soil = np.array([[True]])  # soil land type mask

# Vegetation / snow (static for this demo; dynamic in full model)
tlai = np.array([[0.0]])  # leaf area index [-], bare ground
tsai = np.array([[0.0]])  # stem area index [-]
frac_sno = np.array([[0.0]])  # snow cover fraction [-], no snow

# Soil hydrology (static for this demo)
watsat = np.array([[0.40]])  # saturated volumetric soil water [m³/m³]
h2osoi_ice = np.array([[0.0]])  # frozen soil water [kg/m²], no ice

# ---------------------------------------------------------------------------
# 4. Pre-compute and print diagnostics (mirrors examples/numpy_box_model.py exactly)
# ---------------------------------------------------------------------------
# Representative demo forcing — same values used throughout the time loop
fv_demo = np.array([[0.60]])  # [m/s]   friction velocity
rho_demo = np.array([[1.2]])  # [kg/m³] air density
h2osoi_vol_demo = np.array([[0.05]])  # [m³/m³] volumetric soil water (dry)
h2osoi_liq_demo = np.array([[10.0]])  # [kg/m²] liquid soil water
obu_demo = np.array([[-100.0]])  # [m]     Obukhov length (unstable BL)

fm_d, gwc_sfc_d, gwc_thr_d = moisture_threshold_factor(h2osoi_vol_demo, watsat, fclay)
u_ft0_d = sl2000_u_ft0_dry(rho_demo)[0, 0]
Feff_d, _, _ = drag_partition(z0a, tlai + tsai, frac_rock, frac_veg)
u_s_d = fv_demo[0, 0] * Feff_d[0, 0]
u_it_d = BIT * u_ft0_d
u_st_d = u_ft0_d * fm_d[0, 0] * np.sqrt(rho_demo[0, 0] / RHO_AIR0)
Cd_d = CD0 * np.exp(-CE * (u_st_d - U_ST0) / U_ST0)
k_d = min(CK * (u_st_d - U_ST0) / U_ST0, KAPPA_MAX)
sigma_d = turbulent_sigma_us(np.array([[u_s_d]]), obu_demo)[0, 0]  # wind-speed std at 0.1 m
eta_d = intermittency_factor(
    np.array([[u_s_d]]),
    np.array([[sigma_d]]),
    np.array([[u_it_d]]),
    np.array([[u_ft0_d * fm_d[0, 0]]]),
)[0, 0]

print("--- MUSICA dust-emission box model (Leung et al. 2023/2024) ---")
print("  Precomputed constants:")
print(
    f"    u*ft0_dry at rho_air0 = {U_FT0_DRY_REF:.4f} m/s  "
    f"(Shao & Lu 2000, Dp = {DP_MED * 1e6:.0f} µm)"
)
print(f"    u*it  (Bit * u*ft0)   = {BIT * U_FT0_DRY_REF:.4f} m/s")
print()
print("  Source-to-sink overlap matrix  ovr_src_snk_mss [3 modes × 4 bins]:")
print("    Bin edges [µm]:  0.1-1.0  |  1.0-2.5  |  2.5-5.0  |  5.0-10.0")
for m in range(DST_SRC_NBR):
    row = OVR_SRC_SNK_MSS[m, :]
    print(
        f"    Mode {m + 1} (mss_frc={MSS_FRC_SRC[m]:.3f}): "
        f"{row[0]:.5f}   {row[1]:.5f}   {row[2]:.5f}   {row[3]:.5f}"
    )
print(
    "    Column sums (bin fractions): "
    + "   ".join(f"{OVR_SRC_SNK_MSS[:, n].sum():.5f}" for n in range(NDST))
)
print()
print(
    f"  Cell diagnostics (fv={fv_demo[0, 0]} m/s | rho={rho_demo[0, 0]} kg/m³ | "
    f"fclay={fclay[0, 0]} | h2osoi_vol={h2osoi_vol_demo[0, 0]}):"
)
print(
    f"    gwc_sfc  = {gwc_sfc_d[0, 0]:.4f} kg/kg  |  "
    f"gwc_thr (0.17+0.14·fclay) = {gwc_thr_d[0, 0]:.4f} kg/kg  |  fm = {fm_d[0, 0]:.4f}"
)
print(f"    fclay' (0.1+0.5·min(fclay,0.2)) = {clay_mass_fraction_leung(fclay)[0, 0]:.4f}")
print(
    f"    u*ft0_dry = {u_ft0_d:.4f} m/s  |  u*it = {u_it_d:.4f} m/s  |  "
    f"u*ft_wet = {u_ft0_d * fm_d[0, 0]:.4f} m/s"
)
print(f"    Feff = {Feff_d[0, 0]:.4f}  |  u*s = {u_s_d:.4f} m/s  |  u*st = {u_st_d:.4f} m/s")
print(f"    U(0.1 m) = {u_s_d * _log_law_factor():.3f} m/s  |  sigma_U(0.1 m) = {sigma_d:.4f} m/s")
print(f"    Cd = {Cd_d:.3e}  |  kappa = {k_d:.3f}  |  eta = {eta_d:.6f}")
print()


# ---------------------------------------------------------------------------
# 5. Helper: surface flux -> MICM emission rate
#    F_p [kg m-2 s-1] / dz [m]        = e_mass [kg m-3 s-1]
#    e_mass / DUST_BIN_MW [kg mol-1]  = e      [mol m-3 s-1]  (MICM source)
# ---------------------------------------------------------------------------
def compute_emission_rates(fv, forc_rho, h2osoi_vol, h2osoi_liq, obu, dz):
    """
    Call quacs.dust.leung.dust_emission for the current time step and
    convert surface fluxes [kg m-2 s-1] to volumetric molar emission rates
    [mol m-3 s-1] for the MICM EMIS.* parameters.

    Parameters
    ----------
    fv          : np.ndarray [x, y]   friction velocity [m/s]
    forc_rho    : np.ndarray [x, y]   near-surface air density [kg/m³]
    h2osoi_vol  : np.ndarray [x, y]   volumetric soil moisture [m³/m³]
    h2osoi_liq  : np.ndarray [x, y]   liquid soil water [kg/m²]
    obu         : np.ndarray [x, y]   Obukhov length [m]
    dz          : float               lowest model layer thickness [m]

    Returns
    -------
    dict  {"EMIS.dust_emis_bin_i": [e_i]}  with e_i in [mol m-3 s-1]
    """
    F_p = dust_emission(
        fv=fv,
        forc_rho=forc_rho,
        h2osoi_vol=h2osoi_vol,
        h2osoi_liq=h2osoi_liq,
        h2osoi_ice=h2osoi_ice,
        watsat=watsat,
        tlai=tlai,
        tsai=tsai,
        frac_sno=frac_sno,
        obu=obu,
        fclay=fclay,
        z0a=z0a,
        frac_rock=frac_rock,
        frac_veg=frac_veg,
        Kc_map=Kc_map,
        is_soil=is_soil,
        ovr_src_snk_mss=OVR_SRC_SNK_MSS,
    )  # [x, y, n_bins], kg m-2 s-1

    return {
        f"EMIS.dust_emis_bin_{i + 1}": [F_p[0, 0, i] / dz / DUST_BIN_MW]  # mol m-3 s-1
        for i in range(n_bins)
    }


# ---------------------------------------------------------------------------
# 6. Build solver and initial state
# ---------------------------------------------------------------------------
solver = musica.MICM(
    mechanism=mechanism,
    solver_type=musica.SolverType.rosenbrock_standard_order,
)

state = solver.create_state(number_of_grid_cells=1)

state.set_conditions([300.0], [101325.0])  # T [K], p [Pa]

state.set_concentrations({sp.name: [0.0] for sp in dust_species})

# ---------------------------------------------------------------------------
# 7. Synthetic dynamic forcing (constant over all steps; mirrors examples/numpy_box_model.py)
# ---------------------------------------------------------------------------
dt = 120.0  # s, model time step
dz = 50.0  # m, lowest-layer thickness (MPAS-A / CLM5)
n_steps = 30  # 30 * 120 s = 1 hour

fv_val = np.array([[0.60]])  # m/s, friction velocity (synthetic, constant)
rho_val = np.array([[1.2]])  # kg/m³, near-surface air density
h2osoi_vol = np.array([[0.05]])  # m³/m³, volumetric soil water (dry)
h2osoi_liq = np.array([[10.0]])  # kg/m², liquid soil water
obu_val = np.array([[-100.0]])  # m, Obukhov length (unstable, convective BL)

# ---------------------------------------------------------------------------
# 8. Pre-compute emission rates (constant forcing → constant rates)
# ---------------------------------------------------------------------------
rates = compute_emission_rates(fv_val, rho_val, h2osoi_vol, h2osoi_liq, obu_val, dz)
state.set_user_defined_rate_parameters(rates)

e_per_bin_mol = np.array(
    [rates[f"EMIS.dust_emis_bin_{i + 1}"][0] for i in range(n_bins)]
)  # [mol m-3 s-1]
e_per_bin = e_per_bin_mol * DUST_BIN_MW  # [kg m-3 s-1]

bin_labels = ["0.1-1 µm", "1-2.5 µm", "2.5-5 µm", "5-10 µm"]
print(
    f"  Dynamic forcing: fv={fv_val[0, 0]:.2f} m/s, "
    f"h2osoi_vol={h2osoi_vol[0, 0]:.2f} m³/m³, "
    f"obu={obu_val[0, 0]:.0f} m"
)
print(f"  Integration: dz={dz:.1f} m,  dt={dt:.0f} s,  n_steps={n_steps}")
print()
print(f"  Total surface flux (4 bins): {e_per_bin.sum() * dz:.4e} kg m-2 s-1")
print()
print("  Per-bin emission rates [kg m-3 s-1]:")
for i in range(n_bins):
    print(f"    Bin {i + 1} ({bin_labels[i]}): {e_per_bin[i]:.4e}")
print(
    f"    Total:              {e_per_bin.sum():.4e}  "
    f"(= {e_per_bin_mol.sum():.4e} mol m-3 s-1 with MW = {DUST_BIN_MW} kg/mol)\n"
)

# ---------------------------------------------------------------------------
# 9. Time loop: integrate with Rosenbrock solver
# ---------------------------------------------------------------------------
times = np.zeros(n_steps + 1)
concs = np.zeros((n_steps + 1, n_bins))  # [kg m-3] (converted from mol m-3)

for step in range(n_steps):
    solver.solve(state, dt)

    state_conc = state.get_concentrations()

    for i in range(n_bins):
        concs[step + 1, i] = state_conc[f"dust_bin_{i + 1}"][0] * DUST_BIN_MW  # mol m-3 -> kg m-3
    times[step + 1] = (step + 1) * dt

# ---------------------------------------------------------------------------
# 10. Print evolution + closed-form check at every step
#     For a pure emission system with constant rates, C(t) = e * t exactly.
# ---------------------------------------------------------------------------
print(f"{'time (min)':>10}  {'total dust [kg/m3]':>22}  {'expected [kg/m3]':>22}  match")
print("-" * 70)

all_pass = True

for step in range(n_steps + 1):
    total_actual = concs[step].sum()
    total_expected = e_per_bin.sum() * times[step]  # [kg m-3]

    ok = np.isclose(total_actual, total_expected, rtol=1e-10)
    all_pass = all_pass and ok

    print(
        f"{times[step] / 60:>10.2f}  {total_actual:>22.6e}  "
        f"{total_expected:>22.6e}  {'ok' if ok else 'FAIL'}"
    )

print()
if all_pass:
    print("All steps match the closed-form solution (rtol=1e-10). PASS")
else:
    print("One or more steps deviated from the closed-form solution. FAIL")
