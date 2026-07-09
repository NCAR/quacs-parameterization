# Plume Rise

`quacs.plumerise` is a small Python reference workflow for converting one
wildfire column into a normalized vertical `layer_fraction[k]` profile. It is
written as a direct-input QUACS contract: the public driver accepts
meteorological arrays, vegetation class, and fire-size variables, then returns
only the layer fraction.

![Plume-rise workflow](diagrams/plumerise_workflow.drawio.png)

## Outer Driver And Steps

The public outer driver is
`compute_wildfire_profile_fraction_driver`. It accepts direct meteorological,
vegetation, and fire-size variables for one wildfire column, calls the three
internal workflow steps, and returns only the normalized vertical
`layer_fraction[k]` profile.

```python
from quacs.plumerise import compute_wildfire_profile_fraction_driver

layer_fraction = compute_wildfire_profile_fraction_driver(
    z,
    p,
    t,
    u,
    v,
    qv,
    vegetation_class,
    fire_size_mean,
    fire_size_std=0.0,
)
```

The outer driver calls:

| Step | Module | Main role | Main output | Unit |
| --- | --- | --- | --- | --- |
| Outer driver | `wildfire_profile_fraction_driver.py` | Pass direct input variables through the workflow and expose the public API. | `layer_fraction[k]` | unitless |
| Step 1 | `pyrocb_flag_driver.py` | Diagnose whether the column follows the pyroCb branch using PFT and firepower proxy calculations. | `pyrocb_flag` | unitless boolean |
| Step 2 | `prm_height_driver.py` | Resolve the baseline PRM height range and optional pyroCb-adjusted height range. | `injectH_base_m`, `injectH_top_m`, `injectH_pyroCb_base_m`, `injectH_pyroCb_top_m` | m AGL |
| Step 3 | `layer_fraction_driver.py` | Convert mass split and injection-height ranges into the final vertical profile. | `layer_fraction[k]` | unitless |

## Inputs And Output

The public driver accepts direct variables, not pre-built column, fire, or
vegetation dictionaries. Input variable names do not include unit suffixes;
units are defined by the contract below.

### Meteorological Profile Inputs

| Field | Required? | Unit | Meaning |
| --- | --- | --- | --- |
| `z` | yes | m AGL | Vertical coordinate for the one-column profile. |
| `p` | yes | hPa | Pressure profile on `z`; used by the PFT calculation. |
| `t` | yes | K | Temperature profile on `z`. |
| `u` | yes | m s^-1 | Zonal wind profile on `z`. |
| `v` | yes | m s^-1 | Meridional wind profile on `z`. |
| `qv` | yes | kg kg^-1 | Water-vapor specific humidity profile on `z`. |

### Fire And Vegetation Inputs

| Field | Required? | Unit | Meaning |
| --- | --- | --- | --- |
| `vegetation_class` | yes | unitless | Vegetation lookup key for the internal mass-split and heat-flux tables. |
| `fire_size_mean` | yes | m^2 | Mean fire area or fire-size proxy for the column. |
| `fire_size_std` | optional | m^2 | Fire-size uncertainty used for the high-bound firepower estimate; default is `0.0`. |

### Public Output

| Field | Unit | Meaning |
| --- | --- | --- |
| `layer_fraction[k]` | unitless | Normalized vertical profile fraction for each model layer; the profile sums to 1 over `k`. |

## Repository Layout

```text
quacs/plumerise/                 Python package and three step drivers
quacs/plumerise/examples/        Notebook example
quacs/plumerise/tests/           Contract tests
quacs/plumerise/diagrams/        Workflow diagram source and PNG
```
