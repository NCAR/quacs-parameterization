from pathlib import Path

import numpy as np
import pytest

from quacs.biogenic.megan import MeganSettings, compute_emissions
from quacs.biogenic.megan.src import megan_model

INPUTS = Path(__file__).resolve().parents[1] / "inputs"


@pytest.fixture
def case():
    species = megan_model.load_species_parameters(INPUTS / "EF_LDF.csv")
    pfts = megan_model.load_pft_parameters(INPUTS / "PFT_Fraction.csv")
    return dict(
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


def test_flux_shape_and_sign(case):
    flux = compute_emissions(**case)
    assert flux.shape == (MeganSettings().n_class,)
    assert np.all(np.isfinite(flux))
    assert np.all(flux >= 0.0)
    assert flux[0] > 0.0  # isoprene


def test_no_light_no_isoprene(case):
    """Isoprene has a light-dependent fraction of 1, so it stops in the dark."""
    case["ppfd"] = 0.0
    case["hour"] = 0.0
    flux = compute_emissions(**case)
    assert flux[0] == pytest.approx(0.0, abs=1e-12)


def test_warmer_leaf_emits_more_isoprene(case):
    cool = compute_emissions(**case)[0]
    case["temperature_k"] += 5.0
    warm = compute_emissions(**case)[0]
    assert warm > cool


def test_bad_fraction_shape_raises(case):
    case["pft_fraction"] = np.ones(3)
    with pytest.raises(ValueError):
        compute_emissions(**case)
