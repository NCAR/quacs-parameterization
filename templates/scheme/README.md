# <Scheme name>

<!-- Copy this file to quacs/<process>/<scheme>/README.md and fill in each field.
     Keep the section headings. The CI test tests/test_structure.py checks them. -->

| | |
|---|---|
| **Status** | implemented |
| **Target class** | <Target parameterization class, from the target list> |
| **Contributors** | <Name (institution)> |
| **Original source** | <Link to the Fortran or other code that this port follows> |

<!-- Status must be one of: implemented, in progress, placeholder. -->

One paragraph that tells what physical or chemical process this scheme represents and what it produces.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | |
| Scheme specifics | |
| External datasets | |
| Information from MPAS-A | |
| Other inputs | |
| Information to MPAS-A | |
| Considerations | |
| Equations and references | |

## Classification

These are the dimensions that the QUACS proposal uses to select target parameterizations.

| Dimension | Value |
|---|---|
| Complexity | simple / complex |
| Solving strategy | discrete state update / supplies rates to the solver |
| Aerosol representation | none / bulk / modal / sectional / any |
| Grid | 0-D / 1-D (column) |
| Interdependence | standalone / needs <other parameterization or host component> |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `temperature` | Air temperature | K | `(ncell,)` | MPAS-A |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `flux` | Emission flux | kg m-2 s-1 | `(ncell,)` | MPAS-A state |

## Usage

```python
from quacs.<process>.<scheme> import compute

result = compute(temperature=...)
```

| Example | What it does |
|---|---|
| `examples/musica_box_model.py` | Runs the scheme in a MUSICA box model |
| `examples/<scheme>.ipynb` | Plots the results |

## Tests

```bash
pytest quacs/<process>/<scheme>
```

## References

- Author, A. (Year). Title. Journal, Volume, Pages. https://doi.org/...
