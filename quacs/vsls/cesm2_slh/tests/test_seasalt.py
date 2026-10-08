import musica
import numpy as np
import pytest

from quacs.vsls.cesm2_slh import (
    SEA_SALT_REACTIONS,
    bromine_depletion_factor,
    build_mechanism,
    mean_molecular_speed,
    rate_parameter_key,
    seasalt_mask,
    seasalt_rate_constants,
)

BASE = dict(
    temperature_k=288.0,
    pressure_pa=1.0e5,
    latitude_deg=40.0,
    calday=180.0,
    sad_seasalt_m2_m3=5.0e-5,
    ocean_ice_fraction=1.0,
)

# CAM mo_usrrxt.F90 hard-coded speed prefactors, 100*sqrt(8R/(pi M)) [cm/s/K^0.5]
CAM_SPEED_PREFACTOR = {"het_ss_0": 1.22e3, "het_ss_1": 1.29e3, "het_ss_2": 1.47e3}


@pytest.mark.parametrize("rxn", SEA_SALT_REACTIONS, ids=lambda r: r.name)
def test_speed_matches_cam_prefactor(rxn):
    prefactor_cm = 100.0 * mean_molecular_speed(1.0, rxn.molar_mass_kg_mol)
    assert prefactor_cm == pytest.approx(CAM_SPEED_PREFACTOR[rxn.name], rel=1e-2)


def test_depletion_factor_north_of_30s():
    assert bromine_depletion_factor(0.0, 1.0) == pytest.approx(0.5)
    assert bromine_depletion_factor(-29.9, 200.0) == pytest.approx(0.5)


def test_depletion_factor_southern_extremes():
    # sin term = -1 at calday 0 -> DF_MAX; sin term = +1 at calday 182.5 -> DF_MIN
    assert bromine_depletion_factor(-45.0, 0.0) == pytest.approx(0.9)
    assert bromine_depletion_factor(-45.0, 182.5) == pytest.approx(0.3)
    df = bromine_depletion_factor(-45.0, np.linspace(0, 365, 50))
    assert np.all((df >= 0.3 - 1e-12) & (df <= 0.9 + 1e-12))


def test_mask_cases():
    assert seasalt_mask(40.0, 1.0, 1.0e5) == 1.0  # ocean, lower troposphere
    assert seasalt_mask(40.0, 0.0, 1.0e5) == 0.0  # land north of 60 S
    assert seasalt_mask(-70.0, 0.0, 1.0e5) == 1.0  # land south of 60 S (Antarctica)
    assert seasalt_mask(40.0, 1.0, 2.5e4) == 0.0  # above 300 hPa


def test_hobr_rate_hand_calculation():
    v = np.sqrt(8.0 * 8.314462618 * 288.0 / (np.pi * 96.91e-3))
    expected = 0.25 * 0.0125 * v * 5.0e-5 * 0.5
    assert float(seasalt_rate_constants(**BASE)["het_ss_2"]) == pytest.approx(expected, rel=1e-12)


def test_scaling_factor_is_linear():
    k1 = seasalt_rate_constants(**BASE)["het_ss_2"]
    k3 = seasalt_rate_constants(**BASE, scaling_factor=3.0)["het_ss_2"]
    assert k3 == pytest.approx(3.0 * k1)


def test_vectorized_shape():
    t = np.array([270.0, 280.0, 290.0])
    rates = seasalt_rate_constants(**{**BASE, "temperature_k": t})
    for k in rates.values():
        assert k.shape == (3,)
        assert np.all(np.diff(k) > 0)  # k increases with T through v_mean


def test_negative_sad_rejected():
    with pytest.raises(ValueError):
        seasalt_rate_constants(**{**BASE, "sad_seasalt_m2_m3": -1.0})


def test_micm_matches_analytic_and_br_budget():
    rates = {n: float(k) for n, k in seasalt_rate_constants(**BASE).items()}
    mechanism, _ = build_mechanism()
    solver = musica.MICM(
        mechanism=mechanism, solver_type=musica.SolverType.rosenbrock_standard_order
    )
    state = solver.create_state(1)
    state.set_conditions(temperatures=[288.0], pressures=[1.0e5])
    c0 = {"HOBr": 2e-10, "BrONO2": 1e-10, "BrNO2": 5e-11, "Br2": 0.0, "BrCl": 0.0}
    state.set_concentrations({s: [v] for s, v in c0.items()})
    state.set_user_defined_rate_parameters({rate_parameter_key(n): [k] for n, k in rates.items()})

    dt, n = 300.0, 288  # 24 h
    for _ in range(n):
        solver.solve(state, dt)
    c = {s: v[0] for s, v in state.get_concentrations().items()}

    consumed = 0.0
    for rxn in SEA_SALT_REACTIONS:
        exact = c0[rxn.reactant] * np.exp(-rates[rxn.name] * n * dt)
        assert c[rxn.reactant] == pytest.approx(exact, rel=1e-4)
        consumed += c0[rxn.reactant] - exact

    assert c["Br2"] == pytest.approx(0.65 * consumed, rel=1e-3)
    assert c["BrCl"] == pytest.approx(0.35 * consumed, rel=1e-3)

    def br(x):
        return x["HOBr"] + x["BrONO2"] + x["BrNO2"] + 2.0 * x["Br2"] + x["BrCl"]

    # Gas-phase Br is NOT conserved: +0.65 Br per reactant from sea-salt bromide
    assert br(c) - br(c0) == pytest.approx(0.65 * consumed, rel=1e-3)
