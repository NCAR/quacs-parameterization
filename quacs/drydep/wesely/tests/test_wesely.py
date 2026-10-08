"""Smoke test for the Wesely (1989) scheme.

This test only checks that the scheme runs. The scientist has not written unit tests yet.
"""

import numpy as np

from quacs.drydep.wesely import run_box_model


def test_run_box_model_runs():
    met = dict(
        sfc_temp=298.0,
        air_temp=296.0,
        pressure_sfc=101325.0,
        pressure_10m=100000.0,
        wind_speed=5.0,
        spec_hum=0.01,
        solar_flux=400.0,
        month=7,
        lat=40.0,
    )
    frac_landuse = np.array([0.05, 0.30, 0.10, 0.25, 0.10, 0.05, 0.05, 0.0, 0.05, 0.05, 0.0])
    run_box_model(["O3", "SO2"], met, frac_landuse, box_height_m=1000.0)
