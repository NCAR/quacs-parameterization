# Wildfire Profile Fraction Workflow Contract

This document defines the first QUACS Python reference contract for a one-column
wildfire profile fraction workflow. MPAS-GOCART2G plume-rise code is background
context only; the authoritative deliverable here is the Python contract.

## Core Choice

The core workflow removes the species dimension completely. It returns only:

```python
layer_fraction = compute_wildfire_profile_fraction_driver(
    z,
    p,
    t,
    u,
    v,
    qv,
    vegetation_class,
    fire_size_mean,
)
```

`layer_fraction[k]` is unitless and normalized over model layers. The outer
driver does not special-case missing or tiny fires. If
`fire_size_mean <= min_fire_size`, Step 1 sets `pyrocb_flag=False`, then Step 2
and Step 3 still build the baseline non-pyroCb profile. Diagnostics remain
internal to the step drivers and are not part of the public outer-driver output.

## Main Drivers

| Flow box | Code module | Main input | Main output |
| --- | --- | --- | --- |
| Outer driver | `wildfire_profile_fraction_driver.py` | direct NumPy variables: `z`, `p`, `t`, `u`, `v`, `qv`, `vegetation_class`, `fire_size_mean`, optional `fire_size_std` | `layer_fraction` |
| Step 1: diagnose pyroCb | `pyrocb_flag_driver.py` | `z`, `p`, `t`, `u`, `v`, `qv`, vegetation class, fire size/std | `pyrocb_flag` |
| Step 2: resolve injection heights | `prm_height_driver.py` | Freitas-ready direct variables: `z`, `p`, `t`, `u`, `v`, `qv`, vegetation class, fire size/std, `pyrocb_flag` | `injectH_base_m`, `injectH_top_m`, `injectH_pyroCb_base_m`, `injectH_pyroCb_top_m` |
| Step 3: build layer fraction | `layer_fraction_driver.py` | `z`, vegetation class, direct injection-height variables, `pyrocb_flag` | `layer_fraction` |

## Data Types

### Public Outer-Driver Input Contract

Input variable names do not carry unit suffixes. Units are part of the input
contract table instead.

| Field | Required? | Unit | Meaning |
| --- | --- | --- | --- |
| `z` | yes | m AGL | Vertical coordinate. Current driver assumes `z` can define layer lower bounds, and appends one final top bound using the last spacing. |
| `p` | yes | hPa | Pressure profile on `z`. Required by the pyrometeopy-style PFT calculation. |
| `t` | yes | K | Column temperature profile. |
| `u` | yes | m s^-1 | Zonal wind profile. |
| `v` | yes | m s^-1 | Meridional wind profile. |
| `qv` | yes | kg kg^-1 | Water-vapor specific humidity profile used by the PFT calculation. |
| `vegetation_class` | yes | unitless | Static vegetation class for internal default table lookup. |
| `fire_size_mean` | yes | m^2 | Mean fire area or fire-size proxy for one profile. |
| `fire_size_std` | optional | m^2 | Fire-size uncertainty used for low/high firepower envelopes. Default is `0.0`. |

`fire_exists` is currently `fire_size_mean > min_fire_size`.

The outer driver does not build `column`, `vegetation`, or `fire` wrapper
objects. It passes direct variables into Step 1, Step 2, and Step 3.

### Step Outputs

Step 1 returns the direct bool variable `pyrocb_flag`. Step 2 returns direct
height variables. Step 3 returns direct `layer_fraction`. The public outer
driver returns only `layer_fraction`.

## Step 1: Diagnose PyroCb

Step 1 contains three local helper boxes in `pyrocb_flag_driver.py`:

| Mini box | Module | Input | Output |
| --- | --- | --- | --- |
| Calculate PFT | `pyrocb_flag_driver.py` | `z`, `p`, `t`, `u`, `v`, `qv` | direct `PFT` |
| Calculate qplume | `pyrocb_flag_driver.py` | `fire_size_mean`, `fire_size_std`, `vegetation_class` | direct high-bound `qplume` |
| Diagnose flag | `pyrocb_flag_driver.py` | `PFT`, `qplume` | `pyrocb_flag` |

Logic:

```text
if fire_size_mean <= min_fire_size:
    pyrocb_flag = False
else:
    PFT = compute_PFT(z, p, t, u, v, qv)
    qplume = compute_qplume(fire_size_mean, fire_size_std, vegetation_class)
    pyrocb_flag = qplume >= PFT
```

The PFT implementation follows the `fire_plumes.pft` structure in
`firelab/pyrometeopy`: entrained mixed layer, free-convection level, pressure
weighted mean wind, and the Tory and Kepert PFT formula. There is no public
Step 1 status output. Invalid meteorological profiles or non-positive PFT
components raise errors from the local PFT helper.

## Step 2: Resolve Injection Heights

Step 2 contains three local helper boxes in `prm_height_driver.py`:

| Mini box | Module | Input | Output |
| --- | --- | --- | --- |
| Idealized PRM | `prm_height_driver.py` | Freitas-ready direct variables: `z`, `p`, `t`, `u`, `v`, `qv`, vegetation class, fire size/std | prescribed `injectH_base_m=3000`, `injectH_top_m=6000` |
| Diagnose tropopause | `prm_height_driver.py` | `z`, `t` | WMO tropopause or fallback height |
| Adjust pyroCb height | `prm_height_driver.py` | `pyrocb_flag`, internally diagnosed tropopause | optional pyroCb base/top |

The workflow does not accept tropopause as a direct input. It diagnoses
tropopause from `z` and `t`. If the WMO lapse-rate diagnostic
cannot find a tropopause, diagnostics record fallback use and the height is
`12000 m AGL`.

For pyroCb cases:

```text
injectH_pyroCb_base_m   = tropopause_height_m - 2000
injectH_pyroCb_top_m    = tropopause_height_m + 2000
```

For non-pyroCb cases, all `injectH_pyroCb_*` fields are `None`.

## Step 3: Build Layer Fraction

Step 3 contains three local helper boxes in `layer_fraction_driver.py`:

| Mini box | Module | Input | Output |
| --- | --- | --- | --- |
| Resolve mass split | `layer_fraction_driver.py` | direct `vegetation_class`, `pyrocb_flag` | direct `smoldering_weight`, `flaming_weight`, `pyrocb_weight` |
| Vertical allocation | `layer_fraction_driver.py` | direct `z`, `injectH_base_m`, `injectH_top_m`, optional `injectH_pyroCb_base_m`, `injectH_pyroCb_top_m` | direct `smoldering_layer_fraction`, `flaming_layer_fraction`, `pyrocb_layer_fraction` |
| Combine profiles | `layer_fraction_driver.py` | direct weights and direct layer-fraction variables | final `layer_fraction[k]` |

Smoldering emissions are assigned entirely to the lowest model layer.
Flaming and pyroCb layer profiles use overlap-thickness normalized top-hat
allocation across internally derived `layer_bounds`.

Non-pyroCb split:

```text
smoldering + flaming = 1.0
```

PyroCb split:

```text
pyroCb = 0.5
remaining 0.5 follows the vegetation split table
```

## Shared Tables

Only two shared tables are shown on the diagram.

### `vegetation_split_table`

The public driver does not accept this table. The current workflow uses the
internal default table.

| Column | Meaning |
| --- | --- |
| key | normalized vegetation class |
| `smoldering_fraction` | table fraction before pyroCb adjustment |
| `flaming_fraction` | table fraction before pyroCb adjustment |
| `prm_veg_group` | compact PRM vegetation group |
| `source_mapping` | provenance or fallback tag |

Every row must sum to 1.0. Unknown vegetation classes use the `default` row.

### `heat_flux_table`

The public driver does not accept this table. The current workflow uses the
internal default table when converting fire size to heat flux.

| Column | Meaning |
| --- | --- |
| key | normalized vegetation class |
| `heat_flux_kw_m2` | heat-flux proxy used to convert fire size into equivalent convective firepower |
| `source_mapping` | provenance or fallback tag |

Step 1 estimates firepower from fire size and the internal heat-flux table.

## Workflow Constants

| Constant | Default | Unit | Meaning |
| --- | ---: | --- | --- |
| `min_fire_size` | `0.0` | m^2 | Fire existence threshold. |
| `injectH_base` | `3000.0` | m AGL | Prescribed baseline PRM base height. |
| `injectH_top` | `6000.0` | m AGL | Prescribed baseline PRM top height. |
| `pyroCb_fraction` | `0.5` | unitless | Fraction assigned to pyroCb layer when `pyrocb_flag=True`. |
| `pyroCb_half_width` | `2000.0` | m | Half width around tropopause for pyroCb redistribution. |
| `fallback_tropopause_height` | `12000.0` | m AGL | Tropopause fallback height. |

Step 1 PFT-specific numerical choices are local to `pyrocb_flag_driver.py`.
Step 2 height-specific numerical choices are local to `prm_height_driver.py`.
Step 3 split-specific numerical choices are local to the mass-split helper.

## Optional Inputs And Simplification Candidates

Current optional outer-driver inputs:

| Optional input | Current role | Can simplify? |
| --- | --- | --- |
| `fire_size_std` | Builds low/high fire-size envelope. | Yes. Set to `0.0` or remove if the first workflow is deterministic. |
| internal default tables | Vegetation split and heat flux are still internal assumptions. | Later decide whether each should become a formal table input. |

Recommended first simplification: consider removing public `fire_size_std` if
uncertainty envelopes are not part of the first workflow. Provider functions
and options are already outside the public outer-driver signature.

## Internal Fallback Notes

The public outer driver has no status output. Internal step helpers may still
keep local fallback notes, such as PFT failure, unknown vegetation fallback, or
tropopause fallback, but those notes do not change the public output contract.
