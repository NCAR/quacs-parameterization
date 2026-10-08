import numpy as np
import pytest

from quacs.tagged_tracers.voc_nox import (
    build_1d_grid,
    build_musica_rate_parameters,
    build_single_cell_tracer_emissions,
    build_total_co_emissions,
    convert_surface_flux_to_concentration_tendency,
    transport_1d_advection_diffusion,
)


def test_tracer_emissions_only_at_source():
    grid = build_1d_grid(30, source_cell=6, receptor_cell=20)
    co_emis = build_total_co_emissions(grid)
    tracer_emis = build_single_cell_tracer_emissions(co_emis, source_cell=6)

    assert np.all(co_emis > 0.0)
    assert tracer_emis[6] == co_emis[6]
    assert np.count_nonzero(tracer_emis) == 1


def test_grid_rejects_cells_outside_the_grid():
    with pytest.raises(ValueError):
        build_1d_grid(10, source_cell=10, receptor_cell=2)
    with pytest.raises(ValueError):
        build_1d_grid(10, source_cell=2, receptor_cell=-1)


def test_rate_parameters_have_musica_names():
    params = build_musica_rate_parameters(np.ones(3), np.zeros(3))
    assert set(params) == {
        "EMIS.CO_emission",
        "EMIS.CO_tracer_emission",
        "LOSS.CO_loss",
        "LOSS.CO_tracer_loss",
    }
    assert np.all(params["LOSS.CO_loss"] == 0.0)


def test_diffusion_conserves_mass_with_zero_flux_boundaries():
    c = np.zeros(20)
    c[10] = 1.0
    c_new = transport_1d_advection_diffusion(c, K_diff=0.1, u_adv=0.0, dx=1.0, dt=1.0)
    assert c_new.sum() == pytest.approx(c.sum())
    assert c_new[9] > 0.0 and c_new[11] > 0.0


def test_advection_moves_mass_downwind():
    c = np.zeros(20)
    c[5] = 1.0
    c_new = transport_1d_advection_diffusion(c, K_diff=0.0, u_adv=0.5, dx=1.0, dt=1.0)
    assert c_new[6] > 0.0
    assert c_new[4] == 0.0


def test_surface_flux_conversion():
    assert convert_surface_flux_to_concentration_tendency(2.0, 4.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        convert_surface_flux_to_concentration_tendency(1.0, 0.0)
