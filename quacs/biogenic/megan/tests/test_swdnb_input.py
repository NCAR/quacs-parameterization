"""SWDNB interface contract tests."""

import inspect

from quacs.biogenic.megan import biogenic_emission_megan_v3 as adapter


def test_public_api_expects_swdnb_not_ppfd():
    parameters = inspect.signature(adapter.biogenic_emission_megan_v3).parameters
    assert "swdnb_w_m2" in parameters
    assert "swdnb_24h_w_m2" in parameters
    assert "ppfd" not in parameters
    assert "p24" not in parameters


def test_canopy_receives_shortwave_directly():
    source = inspect.getsource(adapter._calculate_canopy_state)
    assert "solar_w_m2 = swdnb_w_m2" in source
    assert "ppfd / options.solar_to_ppfd" not in source


def test_historical_ppfd_conversion_once():
    source = inspect.getsource(adapter.biogenic_emission_megan_v3)
    assert "float(swdnb_24h_w_m2) * settings.solar_to_ppfd" in source
