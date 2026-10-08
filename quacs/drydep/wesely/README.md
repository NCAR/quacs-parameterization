# drydep/wesely: Wesely (1989) gas dry deposition

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dry Deposition |
| **Contributors** | Saeideh Mohammadi (University of Iowa) |
| **Original source** | CAM `src/chemistry/mozart/mo_drydep.F90`, subroutine `drydep_xactive` (https://github.com/ESCOMP/CAM) |

This scheme computes the dry deposition velocity of gases with the Wesely (1989) resistance model.
The deposition velocity is Vd = 1 / (Ra + Rb + Rc).
Ra is the aerodynamic resistance, Rb is the quasi-laminar sublayer resistance, and Rc is the surface resistance.
Rc comes from 4 parallel pathways: stomata, upper canopy cuticle, lower canopy, and ground.
The scheme computes the resistances for each of the 11 Wesely land types and then makes an area-weighted grid value.
The box model converts Vd to a first-order loss rate k = Vd / H.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Removal of pollutants from the atmosphere, air quality |
| Scheme specifics | Wesely (1989) dry deposition scheme, as in CAM `mo_drydep.F90` (`drydep_xactive`). Ra and Rb follow Walcek et al. (1986). |
| External datasets | Land-use fractions for the 11 Wesely land types. The Wesely resistance tables (5 seasons by 11 land types) are in `wesely.py`. |
| Information from MPAS-A | Surface temperature, air temperature, surface pressure, pressure at the lowest level, 10-m wind speed, specific humidity, surface shortwave flux, snow depth, soil moisture, rain rate |
| Other inputs | Species properties: effective Henry's law constant at 298 K, its temperature dependence, reactivity factor f0, molar mass (`SPECIES_TABLE`). Month and latitude to select the season. Box height H. |
| Information to MPAS-A | Dry deposition velocity for each species. First-order loss rate k = Vd / H for each species. |
| Considerations | The target list asks if Vd belongs in "physics" or "chemistry". In GOCART-2G, the PBL scheme does the mixing and deposition, but the chemistry code computes Vd. This implementation selects the season from the month and the latitude. CAM uses LAI data for this. See the open issues below. |
| Equations and references | Vd = 1 / (Ra + Rb + Rc); Rc = 1 / (1/Rsmx + 1/Rlu + 1/(Rdc + Rcl) + 1/(Rac + Rgs)). Wesely (1989); Walcek et al. (1986). |

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
These are the arguments of `wesely_gas`. `run_box_model` takes the same meteorological values in its `met` dictionary.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `species_name` | Species name, a key of `SPECIES_TABLE` | - | scalar | constant |
| `sfc_temp` | Surface temperature | K | scalar | MPAS-A |
| `air_temp` | Air temperature at the lowest level (about 10 m) | K | scalar | MPAS-A |
| `pressure_sfc` | Surface pressure | Pa | scalar | MPAS-A |
| `pressure_10m` | Pressure at the lowest level (about 10 m) | Pa | scalar | MPAS-A |
| `wind_speed` | 10-m wind speed | m s-1 | scalar | MPAS-A |
| `spec_hum` | Specific humidity | kg kg-1 | scalar | MPAS-A |
| `solar_flux` | Shortwave flux at the surface | W m-2 | scalar | MPAS-A |
| `frac_landuse` | Area fraction of each of the 11 Wesely land types (`LAND_NAMES`) | 1 | `(11,)` | external dataset |
| `month` | Month, 1 to 12 | - | scalar | MPAS-A |
| `lat` | Latitude | degrees | scalar | MPAS-A |
| `snow` | Snow depth (default 0) | m | scalar | MPAS-A |
| `soilw` | Soil moisture fraction (default 0.2) | 1 | scalar | MPAS-A |
| `rain` | Rain rate (default 0) | m s-1 | scalar | MPAS-A |
| `box_height_m` | Box height H for k = Vd / H (default 1000) | m | scalar | MPAS-A (layer thickness) |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the keys of the dictionary that `wesely_gas` returns.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `Vd` | Dry deposition velocity | cm s-1 | scalar | diagnostic |
| `Vd_ms` | Dry deposition velocity | m s-1 | scalar | diagnostic |
| `k` | First-order loss rate, Vd / H | s-1 | scalar | MICM rate parameter |
| `Ra`, `Rb`, `Rc` | Grid-average resistances | s m-1 | scalar | diagnostic |
| `ra_lt`, `rb_lt`, `rc_lt` | Resistances for each land type | s m-1 | `(11,)` | diagnostic |
| `pathways` | Rc pathway resistances for each land type | s m-1 | dictionary of `(11,)` | diagnostic |
| `season`, `heff`, `has_dew`, `has_rain`, `lifetime`, `ustar` | Other diagnostics | various | scalar or `(11,)` | diagnostic |

## Usage

```python
import numpy as np
from quacs.drydep.wesely import wesely_gas

frac = np.zeros(11)
frac[3] = 1.0  # deciduous forest
out = wesely_gas(
    "O3",
    sfc_temp=298.0,
    air_temp=296.0,
    pressure_sfc=101325.0,
    pressure_10m=100000.0,
    wind_speed=5.0,
    spec_hum=0.01,
    solar_flux=400.0,
    frac_landuse=frac,
    month=7,
    lat=40.0,
)
print(out["Vd"], out["k"])
```

| Example | What it does |
|---|---|
| `examples/box_model.py` | Prints the O3 resistances, runs the exponential-decay box model for 8 species, and prints the O3 resistances for each land type |
| `examples/musica_box_model.py` | Gives k = Vd / H to a MUSICA MICM box model and checks the result against the exact solution |
| `examples/wesely_drydep.ipynb` | Plots the concentration decay, Vd and resistances by species, O3 Vd by land type, and the O3 seasonal cycle |

## Tests

```bash
pytest quacs/drydep/wesely
```

## Open issues

- The season comes from the month and the latitude. CAM uses LAI data for each grid cell.
- `HNO3_USE_WESELY_EFFECTIVE_HENRY` is `False`. Thus HNO3 uses the Sander (2015) Henry's law constant, and Rc for HNO3 is near 125 s m-1. Wesely (1989) uses 1e14 M atm-1, which makes HNO3 deposit near the aerodynamic limit.
- The comment on the winter row of `rcls` says that columns 3 and 10 changed from 9000 to 9999. The values in the table at these columns are still 9000. The author must check this row against Wesely (1989), Table 1.
- The functions operate on one grid cell. A vectorized form is necessary for many cells.

## References

- Wesely, M. L. (1989). Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models. Atmos. Environ., 23, 1293-1304. https://doi.org/10.1016/0004-6981(89)90153-4
- Walcek, C. J., Brost, R. A., Chang, J. S., and Wesely, M. L. (1986). SO2, sulfate and HNO3 deposition velocities computed using regional landuse and meteorological data. Atmos. Environ., 20, 949-964. https://doi.org/10.1016/0004-6981(86)90279-9
- Sander, R. (2015). Compilation of Henry's law constants (version 4.0) for water as solvent. Atmos. Chem. Phys., 15, 4399-4981. https://doi.org/10.5194/acp-15-4399-2015
