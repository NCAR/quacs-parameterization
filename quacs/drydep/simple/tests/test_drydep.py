import numpy as np
import pytest

from quacs.drydep.simple.box_model import DEFAULT_MET, compute_drydep_rate
from quacs.drydep.simple.drydep_functions import DIFFG, GET_OBK, METERO


def test_get_obk_sign():
    # Positive heat flux → unstable boundary layer → negative OBK
    obk = GET_OBK(TS=302.0, USTAR=0.15, AIRDEN=1.15, HFLUX=16.8)
    assert np.isfinite(obk)
    assert obk < 0


def test_get_obk_zero_hflux():
    # Zero heat flux must not raise ZeroDivisionError
    obk = GET_OBK(TS=302.0, USTAR=0.15, AIRDEN=1.15, HFLUX=0.0)
    assert np.isfinite(obk)


def test_metero_midpoint_height():
    CZ1, _, _, _ = METERO(130.0, 0.1, 302.0, 0.15, 1.15, 16.8, 0.1, -0.05)
    assert CZ1 == pytest.approx(65.0)


def test_metero_snow_detection():
    _, LSNOW_on, _, _ = METERO(130.0, 0.5, 302.0, 0.15, 1.15, 16.8, 0.1, -0.05)
    _, LSNOW_off, _, _ = METERO(130.0, 0.1, 302.0, 0.15, 1.15, 16.8, 0.1, -0.05)
    assert LSNOW_on
    assert not LSNOW_off


def test_metero_wind_speed():
    _, _, _, W10 = METERO(130.0, 0.1, 302.0, 0.15, 1.15, 16.8, 3.0, 4.0)
    assert W10 == pytest.approx(5.0)  # sqrt(3² + 4²)


def test_diffg_positive():
    assert DIFFG(302.0, 1e5, 201e-3) > 0


def test_diffg_heavier_molecule_slower():
    d_light = DIFFG(302.0, 1e5, 28e-3)
    d_heavy = DIFFG(302.0, 1e5, 201e-3)
    assert d_light > d_heavy


def test_compute_drydep_rate_physical_range():
    dvel, k = compute_drydep_rate(box_height_m=1.0)
    assert dvel > 0
    assert k > 0
    # Hg0 over tropical forest: typically 0.01–0.5 cm/s
    assert 0.001 < dvel < 5.0


def test_compute_drydep_rate_height_scaling():
    dvel1, k1 = compute_drydep_rate(box_height_m=1.0)
    dvel2, k2 = compute_drydep_rate(box_height_m=2.0)
    assert dvel1 == pytest.approx(dvel2)
    assert k1 == pytest.approx(2 * k2)


def test_compute_drydep_rate_met_override():
    _, k_default = compute_drydep_rate()
    _, k_cold = compute_drydep_rate(met={"TC0": 270.0})
    assert k_default != pytest.approx(k_cold)


def test_compute_drydep_rate_partial_override():
    _, k1 = compute_drydep_rate(met={"USTAR": 0.30})
    _, k2 = compute_drydep_rate(met={"USTAR": 0.30, **DEFAULT_MET})
    assert k1 == pytest.approx(k2)
