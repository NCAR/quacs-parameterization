"""Smoke test for the MEGAN3 scheme.

This test only checks that the scheme runs. The scientist has not written unit tests yet.
"""

from pathlib import Path

from quacs.biogenic.megan import compute_emissions
from quacs.biogenic.megan.src import megan_model

INPUTS = Path(__file__).resolve().parents[1] / "inputs"


def test_compute_emissions_runs():
    species = megan_model.load_species_parameters(INPUTS / "EF_LDF.csv")
    pfts = megan_model.load_pft_parameters(INPUTS / "PFT_Fraction.csv")
    compute_emissions(
        day=180.0,
        hour=12.0,
        latitude_deg=45.4,
        temperature_k=298.15,
        pressure_pa=101325.0,
        wind_m_s=2.0,
        relative_humidity_percent=60.0,
        ppfd=1000.0,
        lai=3.0,
        previous_lai=2.8,
        pft_fraction=pfts.fractions_percent / 100.0,
        emission_factor=species.emission_factors_nmol_m2_s,
        light_dependent_fraction=species.light_dependent_fraction,
        t24_k=295.15,
        p24=600.0,
        t10d_k=294.15,
        tmax_k=303.15,
        tmin_k=289.15,
        windmax_m_s=4.0,
    )
