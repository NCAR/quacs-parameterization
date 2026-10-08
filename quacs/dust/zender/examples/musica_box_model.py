# -*- coding: utf-8 -*-
"""
musica_box_model.py
===================
MUSICA box-model for the Zender et al. (2003a) / DEAD dust emission scheme
(CLM4 implementation).

The physics lives in quacs/dust/zender/zender.py.

Forcing matched to the Ginoux box model (quacs/dust/ginoux)
--------------------------------------------------------------
Both box models are driven by the same atmospheric state so their total
surface fluxes can be compared fairly:

    Primary forcing:  w10m = 12.0 m/s  (shared with Ginoux)
    Soil moisture:    gwet = 0.10      (shared with Ginoux)

Zender requires friction velocity fv and 10-m wind u10 rather than w10m
directly. These are derived from w10m via the neutral log-law:

    fv  = κ · w10m / ln(z_ref / z0)
        = 0.4 · 12.0 / ln(10 / 1e-4)
        = 0.4169 m/s              (z0 = 1e-4 m, bare desert soil)

    u10 = w10m = 12.0 m/s        (u10 is defined at the same height)

Ginoux uses gravimetric soil moisture gwet [-]; Zender uses volumetric
h2osoi_vol [m³/m³]. The two are related by:

    gwc_sfc = h2osoi_vol · ρ_water / bd      (Zender Step 2)
    bd      = (1 − watsat) · ρ_quartz

Setting gwc_sfc = gwet = 0.10 and solving:

    h2osoi_vol = gwet · bd / ρ_water
               = 0.10 · 1620 / 1000
               = 0.162 m³/m³

With these mappings both schemes receive identical near-surface wind stress
and identical effective soil moisture; any remaining flux difference reflects
genuine differences in scheme physics (White 1979 saltation law, Fécan
moisture correction, sandblasting efficiency, bin partition).

Units
-----
    F_p                [kg m-2 s-1]   surface flux per bin (science module)
    F_p / dz           [kg m-3 s-1]   volumetric mass source
    F_p / dz / MW      [mol m-3 s-1]  what MICM's EMIS rate expects
    MICM concentration [mol m-3]      -> x MW -> [kg m-3] for output

LOSS:  state.set_user_defined_rate_parameters({"LOSS.Hg0_drydep": [k]})
EMIS:  state.set_user_defined_rate_parameters({"EMIS.dust_emis_bin_i": [e_i]})

References:
    Zender, C. S., Bian, H. S., and Newman, D. (2003a). Mineral Dust
    Entrainment and Deposition (DEAD) model. J. Geophys. Res., 108, 4416.
    https://doi.org/10.1029/2002JD002775

    Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric
    dust cycle. J. Geophys. Res., 100, 16415-16430.
    https://doi.org/10.1029/95JD00690
"""

import musica
import musica.mechanism_configuration as mc
import numpy as np

from quacs.dust.zender.zender import (
    CST_SLT,
    FCLAY_MAX,
    FLX_MSS_FDG_FCT,
    GRAV,
    NDST,
    OVR_SRC_SNK_MSS,
    TMP1,
    dust_emission,
    gwc_threshold,
)

# ---------------------------------------------------------------------------
# 1. Scheme metadata
# ---------------------------------------------------------------------------
n_bins = NDST  # 4 transport bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm
print(f"DEAD / Zender scheme: {n_bins} transport bins")

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
    name="quacs_dust_zender",
    species=dust_species,
    phases=[gas],
    reactions=dust_emissions,
)

# ---------------------------------------------------------------------------
# 3. Static inputs (CLM4 analogs, single cell [1, 1])
# ---------------------------------------------------------------------------
fclay = np.array([[0.10]])  # clay mass fraction [-], 10 %  (matches Ginoux)
mbl_bsn_fct = np.array([[1.0]])  # basin erodibility factor [-] (CLM4 default)
is_soil = np.array([[True]])  # soil land type mask

watsat = np.array([[0.40]])  # saturated volumetric soil water [m³/m³]
h2osoi_ice = np.array([[0.0]])  # frozen soil water [kg/m²]
h2osoi_liq = np.array([[10.0]])  # liquid soil water [kg/m²]  → liqfrac ≈ 1

tlai = np.array([[0.0]])  # leaf area index [-]  (bare ground, matches Ginoux)
tsai = np.array([[0.0]])  # stem area index [-]
frac_sno = np.array([[0.0]])  # snow cover fraction [-]

# ---------------------------------------------------------------------------
# 4. Forcing derivation: w10m → fv, u10  (neutral log-law)
#
#    Ginoux primary input:  w10m = 12.0 m/s
#    Log-law:               fv = κ · w10m / ln(z_ref / z0)
#    Parameters:            κ = 0.4, z_ref = 10 m, z0 = 1e-4 m (bare desert)
#    u10 = w10m (same reference height as Ginoux w10m)
#
#    Soil moisture mapping:
#    Ginoux gwet = 0.10  →  h2osoi_vol = gwet · bd / ρ_water = 0.162 m³/m³
#    (ensures gwc_sfc = gwet in Zender's Step 2, identical moisture effect)
# ---------------------------------------------------------------------------
W10M = 12.0  # m/s — primary wind forcing, shared with Ginoux
GWET = 0.10  # [-]  — primary soil moisture, shared with Ginoux
FORC_RHO = 1.2  # kg/m³

# Neutral log-law constants
KAPPA = 0.4  # von Kármán constant
Z_REF = 10.0  # m, reference height
Z0 = 1.0e-4  # m, aeolian roughness for bare desert soil

fv_derived = KAPPA * W10M / np.log(Z_REF / Z0)  # 0.4169 m/s

# Soil moisture: h2osoi_vol such that gwc_sfc == GWET
RHO_WATER = 1000.0
RHO_QUARTZ = 2700.0
bd = (1.0 - float(watsat[0, 0])) * RHO_QUARTZ  # 1620 kg/m³
h2osoi_vol_eq = GWET * bd / RHO_WATER  # 0.162 m³/m³

# Pack into arrays
fv_val = np.array([[fv_derived]])
u10_val = np.array([[W10M]])  # u10 == w10m at the same height
forc_rho_val = np.array([[FORC_RHO]])
h2osoi_vol = np.array([[h2osoi_vol_eq]])


# ---------------------------------------------------------------------------
# 5. Helper: surface flux -> MICM emission rate
#    F_p [kg m-2 s-1] / dz [m]        = e_mass [kg m-3 s-1]
#    e_mass / DUST_BIN_MW [kg mol-1]  = e      [mol m-3 s-1]  (MICM source)
# ---------------------------------------------------------------------------
def compute_emission_rates(fv, u10, forc_rho, h2osoi_vol, dz):
    """
    Call quacs.dust.zender.dust_emission and convert surface fluxes
    [kg m-2 s-1] to volumetric molar emission rates [mol m-3 s-1].
    """
    F_p = dust_emission(
        fv=fv,
        u10=u10,
        forc_rho=forc_rho,
        h2osoi_vol=h2osoi_vol,
        h2osoi_liq=h2osoi_liq,
        h2osoi_ice=h2osoi_ice,
        watsat=watsat,
        tlai=tlai,
        tsai=tsai,
        frac_sno=frac_sno,
        fclay=fclay,
        mbl_bsn_fct=mbl_bsn_fct,
        is_soil=is_soil,
        tmp1=TMP1,
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
# 7. Integration parameters  (identical to Ginoux box model)
# ---------------------------------------------------------------------------
dt = 120.0  # s
dz = 50.0  # m
n_steps = 30  # 30 × 120 s = 1 hour

# ---------------------------------------------------------------------------
# 8. Pre-compute emission rates and print diagnostics
# ---------------------------------------------------------------------------
# Informational: reproduce Fécan factor and threshold (computed inside
# dust_emission) using the same gwc_threshold() as the science module.
gwc_sfc = h2osoi_vol_eq  # = GWET by construction
gwc_thr = float(gwc_threshold(fclay)[0, 0])  # 0.17 + 0.14*fclay (Zender/CLM)
excess = max(gwc_sfc - gwc_thr, 0.0)
fm = np.sqrt(1.0 + 1.21 * (100.0 * excess) ** 0.68) if gwc_sfc > gwc_thr else 1.0
wnd_frc_thr_slt = TMP1 / np.sqrt(FORC_RHO) * fm  # wet threshold u* [m/s]

# Owen's effect: 10-m threshold wind, then boosted u*
wnd_rfr_thr_slt = W10M * wnd_frc_thr_slt / fv_derived
if W10M >= wnd_rfr_thr_slt:
    wnd_frc_slt = fv_derived + 0.003 * (W10M - wnd_rfr_thr_slt) ** 2
else:
    wnd_frc_slt = fv_derived

# White (1979) horizontal flux and MB95 sandblasting (diagnostic only)
if wnd_frc_slt > wnd_frc_thr_slt:
    r = wnd_frc_thr_slt / wnd_frc_slt
    liq = h2osoi_liq[0, 0] / (h2osoi_ice[0, 0] + h2osoi_liq[0, 0] + 1.0e-6)
    Q_d = (
        CST_SLT * FORC_RHO * wnd_frc_slt**3 / GRAV * (1 - r) * (1 + r) ** 2 * FLX_MSS_FDG_FCT * liq
    )  # [kg m-1 s-1]
else:
    Q_d = 0.0
alpha_d = 100.0 * 10.0 ** (13.4 * min(fclay[0, 0], FCLAY_MAX) - 6.0)  # [m-1]

print("--- MUSICA dust-emission box model (DEAD / Zender et al. 2003a) ---")
print(f"  Forcing matched to Ginoux box model (w10m={W10M} m/s, gwet={GWET}):")
print(f"    fv  (log-law, z0={Z0:.0e} m) = {fv_derived:.5f} m/s")
print(f"    u10 (= w10m)                 = {W10M:.1f} m/s")
print(f"    h2osoi_vol (= gwet·bd/ρ_w)  = {h2osoi_vol_eq:.4f} m³/m³")
print(f"    forc_rho                     = {FORC_RHO:.1f} kg/m³")
print("  Threshold diagnostics:")
print(f"    TMP1                         = {TMP1:.6f} kg^0.5 m^-0.5 s^-1")
print(f"    gwc_sfc = gwet               = {gwc_sfc:.4f} kg/kg")
print(f"    gwc_thr (0.17+0.14·fclay)    = {gwc_thr:.4f} kg/kg")
print(f"    fm (Fécan factor)            = {fm:.4f}")
print(f"    u*_thr (wet)                 = {wnd_frc_thr_slt:.5f} m/s")
print(f"    u10_thr (Owen's ref)         = {wnd_rfr_thr_slt:.4f} m/s")
print(f"    u*_slt (Owen's boosted)      = {wnd_frc_slt:.5f} m/s")
print(f"    Q (White 1979)               = {Q_d:.4e} kg m-1 s-1")
print(f"    alpha (MB95 sandblasting)    = {alpha_d:.4e} m-1")
print(f"    F_tot = Q·alpha              = {Q_d * alpha_d:.4e} kg m-2 s-1")
print(f"  dz={dz:.0f} m, dt={dt:.0f} s, n_steps={n_steps}\n")

rates = compute_emission_rates(fv_val, u10_val, forc_rho_val, h2osoi_vol, dz)
state.set_user_defined_rate_parameters(rates)

e_per_bin_mol = np.array(
    [rates[f"EMIS.dust_emis_bin_{i + 1}"][0] for i in range(n_bins)]
)  # [mol m-3 s-1]
e_per_bin = e_per_bin_mol * DUST_BIN_MW  # [kg m-3 s-1]

bin_labels = ["0.1-1 µm", "1-2.5 µm", "2.5-5 µm", "5-10 µm"]
print(f"Total surface flux (4 bins): {e_per_bin.sum() * dz:.4e} kg m-2 s-1\n")
print("Per-bin emission rates [kg m-3 s-1]:")
for i in range(n_bins):
    print(f"  Bin {i + 1} ({bin_labels[i]}): {e_per_bin[i]:.4e}")
print(
    f"  Total:              {e_per_bin.sum():.4e}  "
    f"(= {e_per_bin_mol.sum():.4e} mol m-3 s-1 with MW = {DUST_BIN_MW} kg/mol)\n"
)

# ---------------------------------------------------------------------------
# 9. Time loop
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
# 10. Print evolution + closed-form check
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
