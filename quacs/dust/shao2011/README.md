# UoC dust emission (Shao et al., 2011, WRF-Chem dust_opt=4)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | WRF-Chem `chem/module_uoc_dust.F` (`uoc_dust_driver`), `chem/module_qf03.F` (`qf03_driver`) and `chem/module_soilpsd.F`, [wrf-model/WRF](https://github.com/wrf-model/WRF), tag v4.8.0, commit `06d4240ae989cc3e50af412bb472df3d9048783c` |

This scheme computes the vertical flux of mineral dust in 5 size bins.
It follows the University of Cologne (UoC) scheme of WRF-Chem (`dust_opt=4`).
Each soil texture class has a particle size distribution for aggregated and for fully dispersed soil.
The friction velocity gives a saltation flux per sand size, and the soil size distribution gives the dust size spectrum.
Erodibility only decides where emission can occur.
The code also keeps the earlier Shao (2001) and Shao (2004) forms (`imod` 1 and 2).

Checked against the compiled pinned Fortran: 166 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms |
| Scheme specifics | Size-resolved dust emission of Shao et al. (2011), simplified form of Shao (2004). The vertical flux is `(1 - cf) * c_y * eta_mi * (1 + sigma_m) * Q * g / u*^2`. The dry threshold follows Shao and Lu (2000). The roughness correction follows Raupach et al. (1993). The moisture correction follows Fécan et al. (1999). The saltation flux follows White (1979). The model option `imod` selects 1 (2001), 2 (2004) or 3 (2011, default). |
| External datasets | Erodibility `erod` (emission where the sum is above 0), area fraction of 12 soil texture classes `soil_frac`. These are static fields. Map them to the MPAS grid. The code holds the soil size distributions. |
| Information from MPAS-A | Friction velocity, volumetric soil moisture, vegetation cover in percent, air density at the surface |
| Other inputs | Size distribution dictionary `psd` from `initialize()`, WRF land surface model number `sf_surface_physics`, option `imod`. Emission coefficient `c_y` and the other constants are in the settings object. |
| Information to MPAS-A | Dust emission flux in 5 bins (upper cut diameters 2, 3.6, 6, 12 and 20 um). Add it to the dust species of the host model. |
| Considerations | Map the 5 bins to the dust species or bins of the host model. The code uses `g = 9.81`. The cover factor `(1 - cf)` acts twice, as in WRF. The code follows WRF v4.8.0, not the version with the fix of Zeng et al. (2020). |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | sectional (5 bins) |
| Grid | 0-D |
| Interdependence | needs friction velocity from the surface layer scheme and soil moisture from the land surface model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `uoc_dust_driver`, the WRF driver of the scheme. It takes the lowest air level and the top soil layer, and takes the soil class fractions from `soil_top_cat`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `dust` | Dust tracers, 5 bins, all levels | ug kg-1 | `(x, level, y, 5)` | MPAS-A |
| `rho_phy` | Air density (the code reads level 1) | kg m-3 | `(x, level, y)` | MPAS-A |
| `dz8w` | Layer thickness (the code reads level 1) | m | `(x, level, y)` | MPAS-A |
| `smois` | Volumetric soil moisture, between 0 and 1 (the code reads soil layer 1) | m3 m-3 | `(x, soil layer, y)` | MPAS-A |
| `ust` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `isltyp` | Soil type index (selects the soil class of `smois_cor`) | 1 | `(x, y)` | external dataset |
| `vegfra` | Vegetation cover | % | `(x, y)` | MPAS-A |
| `dt` | Time step | s | scalar | MPAS-A |
| `emis_dust` | Accumulated emitted dust per bin | kg m-2 | `(x, y, 5)` | MPAS-A |
| `soil_top_cat` | Area fraction of each soil texture class (the code reads the first 12 of the 12 or 16 classes) | 1 | `(x, 12 or 16, y)` | external dataset |
| `erod` | Erodibility (3 parts, only the sum above 0 matters) | 1 | `(x, y, 3)` | external dataset |
| `sf_surface_physics` | WRF land surface model number that selects the dry soil table (default 2, Noah) | 1 | scalar | constant |
| `imod` | Scheme version: 1 Shao 2001, 2 Shao 2004, 3 Shao 2011 (default 3) | 1 | scalar | constant |
| `psd` | Soil particle size distributions, built by `initialize()` from embedded tables. Built inside if omitted. | - | dict | constant |
| `settings` | Tunable source coefficients. See `tunable_names(Shao2011Settings)`. | - | - | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the keys of the dict that `uoc_dust_driver` returns. The driver does not return the flux. It adds flux * `dt` to `emis_dust`, so the flux is (returned `emis_dust` minus the `emis_dust` you passed in) divided by `dt`.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `emis_dust` | Accumulated emitted dust per bin. Divide the change by `dt` to get the flux per bin. | kg m-2 | `(x, y, 5)` | MICM rate parameter |
| `dust` | Dust tracers after the emission (all levels kept) | ug kg-1 | `(x, level, y, 5)` | MPAS-A state |
| `smois` | Soil moisture with the top layer capped at 1 | m3 m-3 | `(x, soil layer, y)` | MPAS-A state |
| `ust_t_min` | Smallest threshold friction velocity over the particle sizes, after the cover and moisture corrections (-9999 where not active) | m s-1 | `(x, y)` | diagnostic |
| `rough_cor` | Roughness correction of the threshold (1 where not active) | 1 | `(x, y)` | diagnostic |
| `smois_cor` | Moisture correction of the threshold, for the soil class `isltyp` (1 where not active) | 1 | `(x, y)` | diagnostic |

The tendency of the lowest layer is dC/dt = flux / (rho_air * dz), with flux = (change of `emis_dust`) / `dt`.

## Usage

```python
import numpy as np
from quacs.dust.shao2011 import initialize, uoc_dust_driver

one = np.ones((1, 1))
level = np.ones((1, 1, 1))
soil_top_cat = np.zeros((1, 12, 1))
soil_top_cat[0, 2, 0] = 1.0
dt = 60.0
out = uoc_dust_driver(
    dust=np.zeros((1, 1, 1, 5)),
    rho_phy=1.1 * level,
    dz8w=10.0 * level,
    smois=0.05 * level,
    ust=0.7 * one,
    isltyp=3 * one.astype(int),
    vegfra=2.0 * one,
    dt=dt,
    emis_dust=np.zeros((1, 1, 5)),
    soil_top_cat=soil_top_cat,
    erod=np.full((1, 1, 3), 0.5),
    psd=initialize(),
)
F_p = out["emis_dust"] / dt  # kg m-2 s-1, shape (1, 1, 5)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the threshold and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 5 bin fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/shao2011
```

## References

- Shao, Y., et al. (2011). Parameterization of size-resolved dust emission and validation with measurements. J. Geophys. Res., 116, D08203. https://doi.org/10.1029/2010JD014527
- Shao, Y. (2004). Simplification of a dust emission scheme and comparison with data. J. Geophys. Res.-Atmos., 109, D10202. https://doi.org/10.1029/2003JD004372
- Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold friction velocity. J. Geophys. Res.-Atmos., 105, 22437-22443. https://doi.org/10.1029/2000JD900304
- White, B. R. (1979). Soil transport by winds on Mars. J. Geophys. Res.-Sol. Ea., 84, 4643-4651. https://doi.org/10.1029/JB084iB09p04643
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157. https://doi.org/10.1007/s00585-999-0149-7
- Raupach, M. R., Gillette, D. A., and Leys, J. F. (1993). The effect of roughness elements on wind erosion threshold. J. Geophys. Res.-Atmos., 98, 3023-3029. https://doi.org/10.1029/92JD01922
- Zeng, Y., et al. (2020). WRF-Chem v3.9 simulations of the East Asian dust storm in May 2017. Geosci. Model Dev., 13, 2125-2147. https://doi.org/10.5194/gmd-13-2125-2020
