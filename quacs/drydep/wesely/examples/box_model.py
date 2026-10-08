#!/usr/bin/env python3
"""
Wesely (1989) gas dry deposition: single-species diagnostic and box model.

This is the demo run from the original notebook by Saeideh Mohammadi. It uses
an exact exponential decay, C(t) = C0 exp(-k t), with k = Vd / H.

Usage::

    python quacs/drydep/wesely/examples/box_model.py
"""

import numpy as np

from quacs.drydep.wesely import LAND_NAMES, N_LANDUSE, run_box_model, wesely_gas

# Mid-latitude summer afternoon
MET = {
    "sfc_temp": 298.0,  # K
    "air_temp": 296.0,  # K
    "pressure_sfc": 101325.0,  # Pa
    "pressure_10m": 100000.0,  # Pa
    "wind_speed": 5.0,  # m/s
    "spec_hum": 0.010,  # kg/kg
    "solar_flux": 400.0,  # W/m2
    "month": 7,
    "lat": 40.0,  # degrees
    "snow": 0.0,  # m
    "soilw": 0.2,  # fraction
    "rain": 0.0,  # m/s
}

# Mixed agricultural and forest landscape, in the order of LAND_NAMES
FRAC_LANDUSE = np.array([0.05, 0.30, 0.10, 0.25, 0.10, 0.05, 0.05, 0.00, 0.05, 0.05, 0.00])

SPECIES = ["O3", "SO2", "NO2", "HNO3", "H2O2", "CO", "CH2O", "NH3"]

if __name__ == "__main__":
    print("\nSINGLE SPECIES DIAGNOSTIC: O3")
    wesely_gas("O3", frac_landuse=FRAC_LANDUSE, verbose=True, **MET)

    print("\nBOX MODEL: ALL SPECIES, 3 HOURS")
    out = run_box_model(
        species_list=SPECIES,
        met=MET,
        frac_landuse=FRAC_LANDUSE,
        box_height_m=1000.0,
        run_hours=3.0,
        time_step=10.0,
    )

    print("\nResistance breakdown by land type for O3:")
    res_o3 = out["results"]["O3"]
    print(
        f"{'Land type':<26} {'frac':<7} {'Ra(s/m)':<10} {'Rb(s/m)':<10} "
        f"{'Rc(s/m)':<10} {'Vd(cm/s)':<10}"
    )
    for lt in range(N_LANDUSE):
        if FRAC_LANDUSE[lt] > 0:
            r_tot = res_o3["ra_lt"][lt] + res_o3["rb_lt"][lt] + res_o3["rc_lt"][lt]
            vd_i = 100.0 / r_tot if r_tot > 0 else 0.0
            print(
                f"{LAND_NAMES[lt]:<26} {FRAC_LANDUSE[lt]:<7.2f} {res_o3['ra_lt'][lt]:<10.1f} "
                f"{res_o3['rb_lt'][lt]:<10.1f} {res_o3['rc_lt'][lt]:<10.1f} {vd_i:<10.4f}"
            )
