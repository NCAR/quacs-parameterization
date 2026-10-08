# drydep/lhs: GEOS-Chem dry deposition with Latin hypercube sampling

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dry Deposition |
| **Contributors** | Ari Feinberg (original Python), QUACS team |
| **Original source** | https://github.com/arifein/offline-drydep |

This scheme is a rewrite of the GEOS-Chem dry deposition code that operates on data, not on indexes. It supports many species, flexible land cover, and many cells at once. A Latin hypercube sampling (LHS) driver explores the sensitivity of the deposition velocity to the inputs for Hg0, SO2, and O3.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Removal of pollutants from the atmosphere, air quality |
| Scheme specifics | GEOS-Chem dry deposition (Wesely, 1989, with GEOS-Chem updates), from the [offline GEOS-Chem dry deposition code](https://github.com/arifein/offline-drydep) |
| External datasets | Olson 2001 land cover and its dry deposition coefficients (`quacs/drydep/data/Olson_2001_Drydep_Inputs.nc`) |
| Information from MPAS-A | Surface temperature, pressure, radiation, cloud fraction, albedo, friction velocity, roughness, 10-m wind, air density, sensible heat flux |
| Other inputs | Species Henry's law constant, reactivity factor, molar mass |
| Information to MPAS-A | Dry deposition velocity and first-order loss rate for each species |
| Considerations | The surface layer is a well-mixed box. The loss rate is k = v_d / H. |
| Equations and references | Wesely (1989); Feinberg et al. (2022) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | none |
| Grid | 0-D |
| Interdependence | standalone |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
The meteorological inputs are keys of the `met` dictionary.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `TC0` | Surface air temperature | K | scalar or `(ncell,)` | MPAS-A |
| `CFRAC` | Cloud fraction | 1 | scalar or `(ncell,)` | MPAS-A |
| `RADIAT` | Incident shortwave radiation | W m-2 | scalar or `(ncell,)` | MPAS-A |
| `AZO` | Roughness height | m | scalar or `(ncell,)` | MPAS-A |
| `USTAR` | Friction velocity | m s-1 | scalar or `(ncell,)` | MPAS-A |
| `PRESSU` | Surface pressure | Pa | scalar or `(ncell,)` | MPAS-A |
| `SUNCOS_MID` | Cosine of the solar zenith angle | 1 | scalar or `(ncell,)` | MPAS-A |
| `ALBD` | Surface albedo | 1 | scalar or `(ncell,)` | MPAS-A |
| `U10M`, `V10M` | 10-m wind components | m s-1 | scalar or `(ncell,)` | MPAS-A |
| `AIRDEN` | Dry air density | kg m-3 | scalar or `(ncell,)` | MPAS-A |
| `HFLUX` | Sensible heat flux | W m-2 | scalar or `(ncell,)` | MPAS-A |
| `box_height_m` | Height of the well-mixed box | m | scalar | MPAS-A (layer thickness) |
| `species` | MUSICA species with `henrys_law_constant` and `reactivity` properties | M atm-1, 1 | list | constant |
| `land_cover` | Olson land-cover patches, each with a fraction and an LAI | 1, m2 m-2 | list, or list of lists `(ncell,)` | external dataset |
| `coefficients` | Olson dry deposition coefficients | various | fixed | external dataset |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `dvel_cms` | Dry deposition velocity (one species, scalar input) | cm s-1 | scalar | diagnostic |
| `k_s` | First-order loss rate (one species, scalar input) | s-1 | scalar | MICM rate parameter |
| `{species: k}` | First-order loss rate for each species (array input) | s-1 | `(ncell,)` for each species | MICM rate parameter |

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

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

Cited in the code:

- Wesely (1989) Atmos. Environ.
- Wang et al. (1998) JGR

From the target parameterization list:

- Wesely, M. L.: Parameterization of Surface Resistances to Gaseous Dry Deposition in Regional-Scale Numerical-Models, Atmos. Environ., 23, 1293–1304, https://doi.org/10.1016/00046981(89)90153-4, 1989.

From the earlier repository README:

- Wesely (1989); Feinberg et al. (2022)
