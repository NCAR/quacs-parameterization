# VOC and NOx tagged tracers (CO prototype)

| | |
|---|---|
| **Status** | in progress |
| **Target class** | VOC + NOx tagged tracers |
| **Contributors** | Liam Sheji (U Albany) |
| **Original source** | New prototype. The target method follows Lapaşcu & Butler (2019). |

This folder holds the first step to the VOC and NOx tagged tracer parameterization.
The prototype tags one species, CO, by one source grid cell.
Total `CO` has emissions in every grid cell, and `CO_tracer` has emissions only in the source cell.
MUSICA solves the local emissions and losses in each cell.
After each MUSICA step, an explicit 1-D advection-diffusion operator moves both species.
The example then shows that `CO_tracer` gets to a receptor cell where its emission is zero.

The next steps are to tag VOC and NOx species, to tag by emission sector and by region, and to use MPAS-A transport.
The 1-D transport operator in this folder is only for the test. It is not a replacement for MPAS-A transport.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, source apportionment |
| Scheme specifics | Gas-phase emissions, photochemistry. This prototype has CO emission and first-order loss only. |
| External datasets | Gridded fire, biogenic, and anthropogenic emissions, initial and boundary conditions. The prototype makes its own emissions with `build_total_co_emissions`. |
| Information from MPAS-A | Chemical species to include, reactions (rates, lifetimes), emission drivers. The prototype uses temperature and pressure. |
| Other inputs | T, RH or q, wind speed and direction, solar radiation, precipitation. The prototype uses a constant wind (`u_adv`) and diffusivity (`K_diff`). |
| Information to MPAS-A | VOC and NOx concentrations tagged by emission sector (anthropogenic, boundary conditions, biomass burning, biogenic, chemistry) or by region (U.S. states, counties). The prototype gives CO and `CO_tracer` tagged by one source cell. |
| Considerations | Boundary condition tracers need spin-up time to mix into the model. Each VOC species must be available to tag and analyze. Input files must have a list of U.S. states and counties with numeric IDs and areas, so that geographic tracers are easy to add. Process analysis must show how transport, mixing, emissions, and chemistry change each tracer. |
| Equations and references | Lapaşcu & Butler (2019) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple (prototype); complex (full VOC and NOx tagging) |
| Solving strategy | supplies rates to the solver (emissions and losses as MICM user-defined rate parameters) |
| Aerosol representation | none |
| Grid | 0-D for the chemistry; the prototype adds 1-D transport for the test |
| Interdependence | needs emissions and host model transport |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `num_grid_cells` | Number of grid cells | 1 | scalar | MPAS-A |
| `source_cell` | Index of the tagged source cell | 1 | scalar | constant |
| `receptor_cell` | Index of the receptor cell for the check | 1 | scalar | constant |
| `temperature` | Air temperature | K | scalar or `(ncell,)` | MPAS-A |
| `pressure` | Air pressure | Pa | scalar or `(ncell,)` | MPAS-A |
| `co_emis` | Total CO emission tendency | mol m-3 s-1 | `(ncell,)` | external dataset (made by `build_total_co_emissions` in the prototype) |
| `co_loss`, `co_tracer_loss` | First-order loss rates of CO and `CO_tracer` (default zero) | s-1 | `(ncell,)` | other parameterization |
| `surface_flux` | Surface emission flux, for conversion to a tendency | mol m-2 s-1 | scalar or `(ncell,)` | external dataset |
| `layer_height` | Layer thickness, for conversion to a tendency | m | scalar or `(ncell,)` | MPAS-A |
| `K_diff` | Diffusivity of the test transport operator | m2 s-1 | scalar | constant |
| `u_adv` | Wind speed of the test transport operator | m s-1 | scalar | constant (MPAS-A in the full scheme) |
| `dx` | Grid spacing | m | scalar | MPAS-A |
| `dt` | Time step | s | scalar | MPAS-A |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `EMIS.CO_emission`, `EMIS.CO_tracer_emission` | Emission rates of CO and `CO_tracer` (from `build_musica_rate_parameters`) | mol m-3 s-1 | `(ncell,)` | MICM rate parameter |
| `LOSS.CO_loss`, `LOSS.CO_tracer_loss` | First-order loss rates of CO and `CO_tracer` | s-1 | `(ncell,)` | MICM rate parameter |
| `co_hist` | CO concentration at the start of each step (from `run_source_receptor`) | mol m-3 | `(ntime, ncell)` | MPAS-A state |
| `co_tracer_hist` | `CO_tracer` concentration at the start of each step | mol m-3 | `(ntime, ncell)` | MPAS-A state |
| proof table | CO and `CO_tracer` at the source and the receptor (from `source_receptor_proof_table`) | mol m-3 | `(ntime,)` columns | diagnostic |

## Usage

```python
import musica
import numpy as np

from quacs.tagged_tracers.voc_nox import (
    build_1d_grid,
    build_co_tracer_mechanism,
    build_musica_rate_parameters,
    build_single_cell_tracer_emissions,
    build_total_co_emissions,
    run_source_receptor,
)

grid = build_1d_grid(30, source_cell=6, receptor_cell=20)
solver = musica.MICM(
    mechanism=build_co_tracer_mechanism(),
    solver_type=musica.SolverType.rosenbrock_standard_order,
)
state = solver.create_state(30)
state.set_conditions(temperatures=np.full(30, 298.0), pressures=np.full(30, 101325.0))
state.set_concentrations({"CO": np.zeros(30), "CO_tracer": np.zeros(30)})

co_emis = build_total_co_emissions(grid)
rates = build_musica_rate_parameters(co_emis, build_single_cell_tracer_emissions(co_emis, 6))
times = np.arange(0, 3010.0, 10.0)
co_hist, co_tracer_hist = run_source_receptor(
    solver, state, rates, times, dt=10.0, K_diff=0.002, u_adv=0.006, dx=1.0
)
```

| Example | What it does |
|---|---|
| `examples/tagged_co_transport.py` | Runs the source-receptor case and checks that `CO_tracer` gets to the receptor |
| `examples/tagged_co_transport.ipynb` | Runs the same case and plots the tracer in space and time |

## Tests

```bash
pytest quacs/tagged_tracers/voc_nox
```

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

Cited in the code:

- None.

From the target parameterization list:

- “Source attribution of European surface O3 using a tagged O3 mechanism” Lapaşcu & Butler, 2019 (https://doi.org/10.5194/acp-19-14535-2019)
