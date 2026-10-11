#!/usr/bin/env python3
"""Single-grid-cell MEGAN3 driver with meteorology defined as Python variables.

No meteorology CSV or plots; EF/LDF output is configured in YAML.
EF/LDF are derived from plant/ecotype/growth-form tables at runtime.
Run with the MEGAN3_in_python repository's ``src`` package importable.

All meteorological numbers below are illustrative, NOT measured site data.
Historical quantities must be supplied by the host model or prepared from
past observations before each call. They are not inferred from one snapshot.
"""
from pathlib import Path
from contextlib import redirect_stdout
import sys
import numpy as np
from ..biogenic_emission_megan_v3 import MeganModel
from ..src.config import MeganConfig

MEGAN_DIR = Path(__file__).resolve().parent.parent

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


def _run_model():
    model = MeganModel.from_yaml(MEGAN_DIR / "configure.yaml")
    Fbio = model.run(
        day=day, hour=hour, latitude_deg=latitude_deg,
        temperature_k=T_k, pressure_pa=P_pa,
        wind_m_s=wind_m_s, relative_humidity_percent=RH_percent,
        swdnb_w_m2=SWDNB, lai=LAI, previous_lai=previous_LAI,
        t24_k=T24_k, swdnb_24h_w_m2=SWDNB_24h, t10d_k=T10d_k,
        tmax_k=Tmax_k, tmin_k=Tmin_k, windmax_m_s=Windmax_m_s,
        kc_7d=kc_7d, swc30d=SWC30D,
    )
    class_names = model.class_names
    n_class = len(class_names)
    np.set_printoptions(precision=6, suppress=True, linewidth=160)
    print(f'MEGAN3: one grid cell, one timestep, {n_class} official MEGAN3.21 classes')
    print('Flux unit: nmol m-2 s-1')
    print(f'day={day:g}, hour={hour:g}')
    for name, flux in zip(class_names, Fbio):
        print(f'  {name}: {flux:.6g}')
    mechanism_flux, species_flux = model.to_mechanism(Fbio, return_species=True)
    print(f'{model.mechanism} fluxes (nmol m-2 s-1):')
    if model.mechanism in ('SPECIES', 'NONE'):
        names = model.mapper.species_names
        for name, flux in zip(names, mechanism_flux):
            print(f'  {name}: {flux:.6g}')
    else:
        for mech_id, name, flux in zip(model.mapper.mechanism_ids, model.mapper.mechanism_names, mechanism_flux):
            print(f'  {mech_id}: {name}: {flux:.6g}')
    print('MGN2MECH mapping diagnostics:', model.mapper.diagnostics())
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

class _Tee:
    """Write stdout to both the terminal and a plain-text result file."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def main():
    config = MeganConfig.from_yaml(MEGAN_DIR / "configure.yaml")
    filename = config.outputs.get('result_file')
    if not isinstance(filename, str) or not filename.strip():
        raise ValueError('files.outputs.result_file must be a nonempty file path')
    output_path = config.path(filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Redirect only during model execution, restoring stdout even if an error occurs.
    with output_path.open('w', encoding='utf-8') as result_stream:
        with redirect_stdout(_Tee(sys.stdout, result_stream)):
            _run_model()


if __name__ == '__main__':
    main()
