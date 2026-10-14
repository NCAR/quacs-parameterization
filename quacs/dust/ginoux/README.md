# GOCART dust emission (Ginoux et al., 2001, WRF-Chem dust_opt=1)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | WRF-Chem `chem/module_gocart_dust.F` (`source_du`, `gocart_dust_driver`) and `phys/module_data_gocart_dust.F`, [wrf-model/WRF](https://github.com/wrf-model/WRF), tag v4.8.0, commit `06d4240ae989cc3e50af412bb472df3d9048783c` |

This scheme computes the vertical flux of mineral dust in 5 size bins.
It follows the GOCART emission of WRF-Chem (`dust_opt=1`).
The flux grows with the 10-m wind speed squared times the excess over a threshold wind speed. A source erodibility map scales it.
The threshold depends on particle size, air density and soil wetness.
Wet soil (saturation fraction 0.5 or more) gets a threshold of 100 m s-1, so it emits nothing.

Checked against the compiled pinned Fortran: 3046 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms |
| Scheme specifics | GOCART dust emission of WRF-Chem v4.8.0 (`dust_opt=1`), after Ginoux et al. (2001). The flux is `C * s_p * S * u10^2 * (u10 - u_t)` for each of 5 bins. The threshold `u_t` depends on size and has a soil-moisture correction. There is no lake term. |
| External datasets | Erodibility `erod` for 3 soil classes (sand, silt, clay). It is a static field. Map it to the MPAS grid. Soil porosity by USGS soil type is a table in the code (used by `compute_gwet`). |
| Information from MPAS-A | 10-m wind speed, soil moisture of the top layer, air density at the surface, land mask |
| Other inputs | Dry threshold wind speed `u_t` per bin from `threshold_velocity()`. Bin radii, bin densities, bin source fractions and the monthly coefficient table `ch_dust` are in the settings object. |
| Information to MPAS-A | Dust emission flux in 5 bins. Add it to the dust species of the host model. |
| Considerations | Map the 5 bins to the dust species or bins of the host model. WRF fixes the month of `ch_dust` at the start of the run. The code uses `g = 9.81`. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | sectional (5 bins) |
| Grid | 0-D |
| Interdependence | standalone (needs `u_t` from `threshold_velocity()` of this scheme and the erodibility map) |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `ginoux_dust_driver`, the WRF driver of the scheme. It takes the lowest air level and the top soil layer.

It sets three inputs of `dust_emission`. The wind speed comes from `u10` and `v10`. If the layer is thinner than 12 m, it comes from `u_phy` and `v_phy`. The threshold wind speed comes from `rho_phy`. The value `gwet` is `smois` over the porosity of the soil type `isltyp`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `dust` | Dust tracers, 5 bins, all levels | ug kg-1 | `(x, level, y, 5)` | MPAS-A |
| `rho_phy` | Air density (the code reads level 1) | kg m-3 | `(x, level, y)` | MPAS-A |
| `dz8w` | Layer thickness (the code reads level 1) | m | `(x, level, y)` | MPAS-A |
| `dt` | Time step | s | scalar | MPAS-A |
| `u10` | 10-m wind, west-east component | m s-1 | `(x, y)` | MPAS-A |
| `v10` | 10-m wind, south-north component | m s-1 | `(x, y)` | MPAS-A |
| `u_phy` | Model wind, west-east component (the code reads level 1 if `dz8w` < 12 m) | m s-1 | `(x, level, y)` | MPAS-A |
| `v_phy` | Model wind, south-north component (the code reads level 1 if `dz8w` < 12 m) | m s-1 | `(x, level, y)` | MPAS-A |
| `smois` | Volumetric soil moisture (the code reads soil layer 1) | m3 m-3 | `(x, soil layer, y)` | MPAS-A |
| `isltyp` | USGS soil type, integer 1 to 19 (selects the porosity) | 1 | `(x, y)` | external dataset |
| `oro` | Land mask (land if `oro < 1.5`) | 1 | `(x, y)` | MPAS-A |
| `erod` | Erodibility for sand, silt, clay | 1 | `(x, y, 3)` | external dataset |
| `emis_dust` | Accumulated emitted dust per bin | kg m-2 | `(x, y, 5)` | MPAS-A |
| `dustin` | Dust tracers handed to other chemistry (used if `chem_opt` is 2 or 11) | ug kg-1 | `(x, y, 5)` | MPAS-A |
| `chem_opt` | WRF chemistry option. 2 and 11 write the update to `dustin`. Default 0. | 1 | scalar | constant |
| `month` | Month that selects the column of `ch_dust` (1 to 12). Default 1. | 1 | scalar | MPAS-A |
| `settings` | Tunable source coefficients. See `tunable_names(GinouxSettings)`. | - | - | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the keys of the dict that `ginoux_dust_driver` returns. Sea cells keep their starting values.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `flux` | Dust emission flux per bin | kg m-2 s-1 | `(x, y, 5)` | MICM rate parameter |
| `dust` | Dust tracers after the emission (all levels kept, unchanged if `chem_opt` is 2 or 11) | ug kg-1 | `(x, level, y, 5)` | MPAS-A state |
| `dustin` | Updated lowest-level tracers if `chem_opt` is 2 or 11 | ug kg-1 | `(x, y, 5)` | MPAS-A state |
| `emis_dust` | Accumulated emitted dust per bin | kg m-2 | `(x, y, 5)` | diagnostic |

The tendency of the lowest layer is dC/dt = flux / (rho_air * dz).

## Usage

```python
import numpy as np
from quacs.dust.ginoux import ginoux_dust_driver

one = np.ones((1, 1))
level = np.ones((1, 1, 1))
out = ginoux_dust_driver(
    dust=np.zeros((1, 1, 1, 5)),
    rho_phy=1.2 * level,
    dz8w=10.0 * level,
    dt=60.0,
    u10=12.0 * one,
    v10=0.0 * one,
    u_phy=12.0 * level,
    v_phy=0.0 * level,
    smois=0.034 * level,
    isltyp=np.array([[1]]),
    oro=one,
    erod=np.full((1, 1, 3), 0.5),
    emis_dust=np.zeros((1, 1, 5)),
    dustin=np.zeros((1, 1, 5)),
)
flux = out["flux"]  # kg m-2 s-1, shape (1, 1, 5)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the bin thresholds and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 5 bin fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/ginoux
```

## References

- Ginoux, P., et al. (2001). Sources and distributions of dust aerosols simulated with the GOCART model. J. Geophys. Res.-Atmos., 106, 20255-20273. https://doi.org/10.1029/2000JD000053
