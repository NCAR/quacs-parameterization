# drydep/lhs

Extends the simple MUSICA dry deposition box model with Latin Hypercube Sampling (LHS) to explore sensitivity across a wide range of atmospheric conditions. Runs Hg0, SO2, and O3 simultaneously and produces sensitivity plots of deposition velocity and loss rate against each input dimension.

This uses the [GEOS-Chem standalone drypdep scheme](https://github.com/arifein/offline-drydep), but
using a refactored version of the code to look more like python and operate on data rather than
on indexes.

## Examples

### `musica_box_model.py` — single-cell reference run

Runs a single tropical mixed-forest grid cell with hardcoded meteorology and land cover. Computes deposition velocity, first-order loss rate, and e-folding lifetime for each species, then integrates concentrations forward with a MUSICA MICM Rosenbrock solver and prints a concentration table.

**Usage**

```
python examples/musica_box_model.py
```

**Configuration** (edit constants at the top of the file)

| Variable | Description |
|----------|-------------|
| `met` | Meteorological inputs (temperature, pressure, radiation, wind, etc.) |
| `land_cover` | List of Olson land-cover patches with fraction and LAI |
| `initial_ng_m3` | Initial concentrations for Hg0, SO2, O3 (ng m⁻³) |

**Output** — printed to stdout: deposition velocities and loss rates, followed by a time series table of concentrations (mol m⁻³) at ~20 time points spanning 3 e-folding times of the fastest-depositing species.

---

### `lhs_driver.py` — LHS ensemble

Generates an ensemble of atmospheric conditions with Latin Hypercube Sampling across meteorological variables, land-cover fractions, land-cover LAI, and initial concentrations. Runs the MUSICA MICM box model for each cell and produces sensitivity plots of deposition velocity against each input dimension.

**Usage**

```
python examples/lhs_driver.py
python examples/lhs_driver.py --cells 200 --seed 42
```

**Arguments**

| Argument | Default | Description |
|----------|---------|-------------|
| `--cells` | 100 | Number of LHS samples |
| `--seed` | 0 | RNG seed for reproducibility |
| `--output` | `output` | Directory for output files |

**Sampled dimensions**

| Group | Dimensions |
|-------|-----------|
| Meteorology | Temperature, pressure, cloud fraction, radiation, albedo, wind speed, friction velocity, sensible heat flux, boundary layer height, air density |
| Land cover | Fraction and LAI for each Olson land-cover archetype |
| Initial concentrations | Hg0 (ng m⁻³), SO2 (ng m⁻³), O3 (ng m⁻³) |

**Outputs**

| File | Description |
|------|-------------|
| `output/lhs_results.csv` | Per-cell inputs, deposition velocities, loss rates, and final concentrations |
| `output/lhs_sensitivity_<species>.png` | Scatter plots of v_d vs. each meteorological input |
| `output/lhs_landcover_<species>.png` | Scatter plots of v_d vs. land-cover fraction |
| `output/lhs_timeseries.png` | Mean ± 1 std normalised concentration over time across all cells |
