# drydep/simple: GEOS-Chem dry deposition for Hg0

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dry Deposition |
| **Contributors** | Ari Feinberg (original Python), QUACS team |
| **Original source** | https://github.com/arifein/offline-drydep |

This scheme is a direct translation of the offline GEOS-Chem dry deposition code for Hg0. It computes a dry deposition velocity from meteorology and Olson land-cover parameters. It converts the velocity to a first-order loss rate. A MUSICA MICM box model then integrates the Hg0 concentration with the Rosenbrock solver.

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
| (land cover) | Tropical mixed forest, fixed in `box_model.py` | 1 | fixed | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `dvel_cms` | Hg0 dry deposition velocity | cm s-1 | scalar | diagnostic |
| `k_s` | Hg0 first-order loss rate | s-1 | scalar | MICM rate parameter |

## Usage

**Simple 0-D calculation** (no MUSICA required):
```
python examples/ex_drydep_simple.py
```
Prints the dry deposition velocity for a single atmospheric column with default tropical mixed-forest conditions.

**MUSICA box model** (requires MUSICA installed):
```
python examples/musica_drydep_box_model.py
```
Integrates Hg0 concentration over ~5 e-folding lifetimes and prints a time series of concentration decay.

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

Cited in the code:

- Wesely [1989]
- Wang et al., JGR, 1998
- Offline GEOS-Chem dry deposition code by Ari Feinberg (`@author: arifeinberg`)

From the target parameterization list:

- Wesely, M. L.: Parameterization of Surface Resistances to Gaseous Dry Deposition in Regional-Scale Numerical-Models, Atmos. Environ., 23, 1293–1304, https://doi.org/10.1016/00046981(89)90153-4, 1989.

From the earlier repository README:

- Wesely (1989); Feinberg et al. (2022)
