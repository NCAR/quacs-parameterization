import numpy as np
import pytest

from quacs.drydep.wesely import SPECIES_TABLE, get_season_index, wesely_gas

MET = dict(
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
FRAC = np.array([0.05, 0.30, 0.10, 0.25, 0.10, 0.05, 0.05, 0.00, 0.05, 0.05, 0.00])


@pytest.mark.parametrize("species", sorted(SPECIES_TABLE))
def test_vd_is_positive_and_finite(species):
    out = wesely_gas(species, frac_landuse=FRAC, **MET)
    assert np.isfinite(out["Vd"])
    assert out["Vd"] > 0.0
    assert out["k"] == pytest.approx(out["Vd_ms"] / 1000.0)


def test_o3_vd_in_plausible_range():
    vd = wesely_gas("O3", frac_landuse=FRAC, **MET)["Vd"]
    assert 0.1 < vd < 2.0  # cm/s, typical daytime summer values


def test_soluble_gas_deposits_faster_than_o3():
    vd_o3 = wesely_gas("O3", frac_landuse=FRAC, **MET)["Vd"]
    vd_so2 = wesely_gas("SO2", frac_landuse=FRAC, **MET)["Vd"]
    assert vd_so2 > vd_o3


def test_night_o3_vd_is_smaller_than_day():
    day = wesely_gas("O3", frac_landuse=FRAC, **MET)["Vd"]
    night = wesely_gas("O3", frac_landuse=FRAC, **{**MET, "solar_flux": 0.0})["Vd"]
    assert night < day


def test_season_index():
    assert get_season_index(7, 40.0) == 0
    assert get_season_index(1, -40.0) == 0  # January is summer in the south
    assert get_season_index(7, 40.0, snow=0.1) == 3


def test_unknown_species_raises():
    with pytest.raises(ValueError):
        wesely_gas("XYZ", frac_landuse=FRAC, **MET)


def test_bad_landuse_shape_raises():
    with pytest.raises(ValueError):
        wesely_gas("O3", frac_landuse=np.ones(5), **MET)
