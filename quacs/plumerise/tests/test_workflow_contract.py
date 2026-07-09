import numpy as np
import pytest

from quacs.plumerise import compute_wildfire_profile_fraction_driver


def synthetic_profile_inputs():
    z = np.arange(0.0, 16000.0, 1000.0)
    t = np.where(z <= 10000.0, 300.0 - 6.5 * (z / 1000.0), 235.0)
    p = 1013.25 * np.exp(-z / 8000.0)
    return {
        "z": z,
        "p": p,
        "t": t,
        "u": np.ones_like(z),
        "v": np.zeros_like(z),
        "qv": np.full_like(z, 0.01),
    }


def no_tropopause_profile_inputs():
    z = np.arange(0.0, 16000.0, 1000.0)
    p = 1013.25 * np.exp(-z / 8000.0)
    return {
        "z": z,
        "p": p,
        "t": 300.0 - 6.5 * (z / 1000.0),
        "u": np.ones_like(z),
        "v": np.zeros_like(z),
        "qv": np.full_like(z, 0.01),
    }


def test_zero_fire_size_still_uses_non_pyrocb_profile():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **synthetic_profile_inputs(),
        vegetation_class="forest",
        fire_size_mean=0.0,
    )

    assert sum(layer_fraction) == pytest.approx(1.0)
    assert layer_fraction[0] == pytest.approx(0.30)
    assert layer_fraction[3] == pytest.approx(0.70 / 3.0)
    assert layer_fraction[4] == pytest.approx(0.70 / 3.0)
    assert layer_fraction[5] == pytest.approx(0.70 / 3.0)


def test_non_pyrocb_uses_default_table_split_and_overlap_allocation():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **synthetic_profile_inputs(),
        vegetation_class="forest",
        fire_size_mean=1000.0,
        fire_size_std=100.0,
    )

    assert sum(layer_fraction) == pytest.approx(1.0)
    assert layer_fraction[0] == pytest.approx(0.30)
    assert layer_fraction[1] == pytest.approx(0.0)
    assert layer_fraction[2] == pytest.approx(0.0)
    assert layer_fraction[3] == pytest.approx(0.70 / 3.0)
    assert layer_fraction[4] == pytest.approx(0.70 / 3.0)
    assert layer_fraction[5] == pytest.approx(0.70 / 3.0)


def test_pyrocb_assigns_half_to_tropopause_range_and_remainder_to_split():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **synthetic_profile_inputs(),
        vegetation_class="forest",
        fire_size_mean=20_000_000.0,
    )

    assert sum(layer_fraction) == pytest.approx(1.0)
    assert layer_fraction[0] == pytest.approx(0.15)
    assert layer_fraction[3] == pytest.approx(0.35 / 3.0)
    assert layer_fraction[4] == pytest.approx(0.35 / 3.0)
    assert layer_fraction[5] == pytest.approx(0.35 / 3.0)
    for idx in range(8, 12):
        assert layer_fraction[idx] == pytest.approx(0.5 / 4.0)


def test_unknown_vegetation_uses_default_tables():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **synthetic_profile_inputs(),
        vegetation_class="unknown_land_cover",
        fire_size_mean=1000.0,
        fire_size_std=100.0,
    )

    assert sum(layer_fraction) == pytest.approx(1.0)
    assert layer_fraction[0] == pytest.approx(0.40)
    assert layer_fraction[3] == pytest.approx(0.60 / 3.0)


def test_pyrocb_tropopause_fallback_uses_default_height():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **no_tropopause_profile_inputs(),
        vegetation_class="forest",
        fire_size_mean=20_000_000.0,
    )

    assert sum(layer_fraction) == pytest.approx(1.0)
    for idx in range(10, 14):
        assert layer_fraction[idx] == pytest.approx(0.5 / 4.0)
