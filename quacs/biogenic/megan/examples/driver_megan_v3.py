#!/usr/bin/env python3
"""Single-grid-cell MEGAN3 driver with meteorology defined as Python variables.

No namelist, meteorology CSV, output files, or plots.
EF/LDF and PFT fractions still use the original upstream CSV readers.
Run with the MEGAN3_in_python repository's ``src`` package importable.

All meteorological numbers below are illustrative, NOT measured site data.
Historical quantities must be supplied by the host model or prepared from
past observations before each call. They are not inferred from one snapshot.
"""
from pathlib import Path
import sys

import numpy as np
from .. import biogenic_emission_megan_v3 as adapter
from ..src.MGN2MECH import MGN2MECH

# ============================================================================
# 1. FILES: only EF/LDF and PFT CSV remain
# ============================================================================
# Directory containing this driver
MEGAN_DIR = Path(__file__).resolve().parent.parent
emission_factors_file = MEGAN_DIR / "inputs" / "EF_LDF.csv"
pft_fractions_file = MEGAN_DIR / "inputs" / "PFT_Fraction.csv"
# ============================================================================
# 2. MODEL / SITE: former namelist settings
# ============================================================================
n_class = 19  # number of VOC classes
NLayers = 5   # number of canopy layers
canopy_layer_method = 'gaussian'  # 'gaussian' or 'uniform'
NRTYP = 6 # number of canopy types

solar_constant_w_m2 = 1361.5
solar_to_ppfd = 2.1
air_quality_index = 40.0
co2_ppm = 403.0
kc_min = 0.0
kc_max = 0.82
wilting_point = 0.196
humidity_mode = 'RH'  # Current emission adapter accepts RH only

# ============================================================================
# 3. ACTIVITY SWITCHES: former namelist settings
# ============================================================================
GAMBD_YN = True
GAMAQ_YN = True
GAMHT_YN = False
GAMLT_YN = False
GAMHW_YN = False
GAMSM_YN = False
GAMCO2_YN = False

# ============================================================================
# 4. METEOROLOGY FOR ONE GRID CELL AND ONE TIMESTEP (SCALARS)
#    Illustrative values only. Update these values before each online call.
# ============================================================================
latitude_deg = 45.4
day = 180.0                            # Day of year
hour = 12.0                            # Hour
T_k = 298.15                           # Air temperature, K at 2m
P_pa = 101325.0                        # Surface air pressure, Pa
wind_m_s = 2.0                         # Wind speed, m s-1 at 2m
RH_percent = 60.0                      # Relative humidity, %
QV_kg_kg = 0.01594                     # water vapor mixing ratio mixing ratio
SWDNB = 500.0                          # MPAS-A downward shortwave radiation at bottom, W m-2
LAI = 3.0                              # Current LAI, m2 m-2
previous_LAI = 2.8                     # Previous LAI (8 days ago)for leaf-age response

# Historical statistics supplied explicitly by the host model.
# QUACS uses causal rolling-24h history; upstream standalone used current-day means.
T24_k = 295.15                         # Original: Rolling previous-24h mean T, K(include current hour)
SWDNB_24h = 200.0                      # Rolling previous-24h mean SWDNB, W m-2,include current hour)
T10d_k = 294.15                        # Mean T of previous 10 complete days, K
Tmax_k = 303.15                        # Original: current-day maximum T, K
Tmin_k = 289.15                        # Original: current-day minimum T, K
Windmax_m_s = 4.0                      # Original: current-day maximum wind, m s-1
SWC30D = 0.25                          # Mean Soil water content of prvious 20 complete days, m3 m-3
kc_7d = 0.2436                         # Mean evapotranspiration ratio of previous 7 complete days, unitless
# ============================================================================
# 5. OPTIONAL MIXING-RATIO UPDATE AND SCREEN PRINT
# ============================================================================
integrate_mixing_ratio = False
bottom_layer_dz_m = 50.0 # depth of lowest model layer
step_seconds = 3600.0 # time step
tracer_initial_mol_mol = 0.0 # initial concentrations of VOC 
ALT_m3_kg = 1.0 / 1.2                # WRF lowest-layer ALT, m3/kg dry air
M_DRY_AIR_KG_MOL = 0.02897  #Average molar mass of dry air

def emission_flux_to_ppmv_tendency(flux_nmol_m2_s, alt_m3_kg, dz_m):
    """Convert surface molar flux to dry-air ppmv/s in the lowest layer.

    Assumes all flux enters and mixes uniformly through the lowest layer;
    does not include transport, chemistry, deposition, or layer mass changes.
    """
    if not np.isfinite(alt_m3_kg) or alt_m3_kg <= 0:
        raise ValueError('ALT_m3_kg must be a positive finite inverse dry-air density.')
    if not np.isfinite(dz_m) or dz_m <= 0:
        raise ValueError('bottom_layer_dz_m must be positive and finite.')
    return np.asarray(flux_nmol_m2_s) * 1.0e-3 * M_DRY_AIR_KG_MOL * alt_m3_kg / dz_m


def main():
    if humidity_mode != 'RH':
        raise ValueError('Current MEGAN adapter supports RH, not specific humidity.')
    if not np.isfinite(SWC30D) or not (0.0 <= SWC30D <= 1.0):
        raise ValueError('SWC30D must be a volumetric soil water content between 0 and 1 m3/m3.')
    if GAMSM_YN and not np.isfinite(kc_7d):
        raise ValueError('Upstream MEGAN GAMSM_YN uses kc_7d, not SWC30D. Provide kc_7d or keep GAMSM_YN=False.')
    mapper = MGN2MECH()
    species = adapter.load_species_parameters(emission_factors_file)
    if tuple(species.names) != tuple(mapper.class_names):
        raise ValueError("EF_LDF.csv classes must match official MEGAN3.21 class order")
    pfts = adapter.load_pft_parameters(pft_fractions_file)
    if len(species.names) != n_class:
        raise ValueError(f'n_class={n_class}, but EF/LDF has {len(species.names)} classes')
    if len(pfts.fractions_percent) != NRTYP:
        raise ValueError(f'NRTYP={NRTYP}, but PFT input has {len(pfts.fractions_percent)} types')

    settings = adapter.MeganSettings(
        n_class=n_class, NLayers=NLayers, NRTYP=NRTYP,
        canopy_layer_method=canopy_layer_method,
        solar_constant_w_m2=solar_constant_w_m2,
        solar_to_ppfd=solar_to_ppfd,
        air_quality_index=air_quality_index,
        co2_ppm=co2_ppm, kc_min=kc_min, kc_max=kc_max,
    )
    switches = dict(
        bidirectional_lai_response=GAMBD_YN,
        air_quality_response=GAMAQ_YN,
        high_temperature_response=GAMHT_YN,
        low_temperature_response=GAMLT_YN,
        high_wind_response=GAMHW_YN,
        soil_moisture_response=GAMSM_YN,
        co2_response=GAMCO2_YN,
    )

    pft_fraction = pfts.fractions_percent / 100.0
    ef = species.emission_factors_nmol_m2_s
    ldf = species.light_dependent_fraction

    Fbio = adapter.biogenic_emission_megan_v3(
        day=day, hour=hour, latitude_deg=latitude_deg,
        temperature_k=T_k, pressure_pa=P_pa,
        wind_m_s=wind_m_s, relative_humidity_percent=RH_percent,
        swdnb_w_m2=SWDNB, lai=LAI, previous_lai=previous_LAI,
        pft_fraction=pft_fraction, emission_factor=ef,
        light_dependent_fraction=ldf,
        t24_k=T24_k, swdnb_24h_w_m2=SWDNB_24h, t10d_k=T10d_k,
        tmax_k=Tmax_k, tmin_k=Tmin_k, windmax_m_s=Windmax_m_s,
        kc_7d=kc_7d, swc30d=SWC30D, settings=settings, switches=switches,
    )
    np.set_printoptions(precision=6, suppress=True, linewidth=160)
    print(f'MEGAN3: one grid cell, one timestep, {n_class} official MEGAN3.21 classes')
    print('Flux unit: nmol m-2 s-1')
    print(f'day={day:g}, hour={hour:g}, SWC30D={SWC30D:g} m3/m3')
    for name, flux in zip(species.names, Fbio):
        print(f'  {name}: {flux:.6g}')
    cbmz_flux, species_flux = mapper.convert(Fbio, 'CBMZ', return_species=True)
    print('CBMZ fluxes (nmol m-2 s-1):')
    for mech_id, name, flux in zip(mapper.mechanism_ids, mapper.mechanism_names, cbmz_flux):
        print(f'  {mech_id}: {name}: {flux:.6g}')
    print('MGN2MECH mapping diagnostics:', mapper.diagnostics())
    if integrate_mixing_ratio:
        if not np.isfinite(step_seconds) or step_seconds < 0:
            raise ValueError('step_seconds must be nonnegative and finite.')
        tendency_ppmv_s = emission_flux_to_ppmv_tendency(
            Fbio, ALT_m3_kg, bottom_layer_dz_m
        )
        x_ppmv = np.full(n_class, tracer_initial_mol_mol * 1.0e6, dtype=float)
        x_ppmv += tendency_ppmv_s * step_seconds
        print('Tracer tendencies (ppmv dry-air/s):', tendency_ppmv_s)
        print('Tracer mixing ratios after one timestep (ppmv dry-air):', x_ppmv)

if __name__ == '__main__':
    main()
