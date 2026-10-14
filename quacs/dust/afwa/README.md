# AFWA dust emission (LeGrand et al., 2019, WRF-Chem dust_opt=3)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | WRF-Chem `chem/module_gocart_dust_afwa.F` (`source_dust`) and `phys/module_data_gocart_dust.F`, [wrf-model/WRF](https://github.com/wrf-model/WRF), tag v4.8.0, commit `06d4240ae989cc3e50af412bb472df3d9048783c` |

This scheme computes the vertical flux of mineral dust in 5 size bins.
It follows the AFWA scheme of WRF-Chem (`dust_opt=3`).
The friction velocity drives a horizontal saltation flux over 9 sand size bins, weighted by the soil texture.
A clay-dependent sandblasting efficiency turns this into a vertical flux.
A fixed lognormal distribution (Kok, 2011) splits the vertical flux into the 5 dust bins.

Checked against the compiled pinned Fortran: 4700 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms |
| Scheme specifics | AFWA dust emission of WRF-Chem v4.8.0 (LeGrand et al., 2019). The dry threshold follows Iversen and White (1982). The moisture correction follows Fécan et al. (1999). The saltation flux follows Marticorena and Bergametti (1995). Emission is zero over water, over soil types 15, 16 and 18, where roughness length is above 0.2 m, and where snow depth is above 0.01 m. |
| External datasets | Erodibility `erod` (3 classes), clay and sand mass fractions, soil type. Optional: vegetation fraction, MODIS vegetation mask, DRI erodibility, NGA soil texture. These are static fields. Map them to the MPAS grid. |
| Information from MPAS-A | Friction velocity, volumetric soil moisture, air density at the surface, roughness length, snow depth, land mask |
| Other inputs | WRF option switches (`dust_veg`, `dust_smois`, `dust_dsr`, `dust_soils`, `sf_surface_physics`). Tuning constants `alpha`, `gamma`, `smtune`, `ustune` and the bin distribution are in the settings object. |
| Information to MPAS-A | Dust emission flux in 5 bins. Add it to the dust species of the host model. |
| Considerations | Map the 5 bins to the dust species or bins of the host model. The code uses `g = 9.81`. The option `sf_surface_physics` must match the land surface model of the host. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | sectional (5 bins) |
| Grid | 0-D |
| Interdependence | needs friction velocity from the surface layer scheme and soil fields from the land surface model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `afwa_dust_driver`, the WRF driver of the scheme. It takes the lowest air level and the top soil layer, and passes the keyword arguments below (`**fields`) to `dust_emission_afwa`. It sets the air density and soil moisture of `dust_emission_afwa` from `rho_phy` and `smois`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `dust` | Dust tracers, 5 bins, all levels | ug kg-1 | `(x, level, y, 5)` | MPAS-A |
| `rho_phy` | Air density (the code reads level 1) | kg m-3 | `(x, level, y)` | MPAS-A |
| `dz8w` | Layer thickness (the code reads level 1) | m | `(x, level, y)` | MPAS-A |
| `smois` | Volumetric soil moisture (the code reads soil layer 1) | m3 m-3 | `(x, soil layer, y)` | MPAS-A |
| `dt` | Time step | s | scalar | MPAS-A |
| `u10` | 10-m wind, west-east component (for the lofting diagnostic) | m s-1 | `(x, y)` | MPAS-A |
| `v10` | 10-m wind, south-north component (for the lofting diagnostic) | m s-1 | `(x, y)` | MPAS-A |
| `emis_dust` | Accumulated emitted dust per bin | kg m-2 | `(x, y, 5)` | MPAS-A |
| `tot_edust` | Accumulated emitted dust, all bins | kg m-2 | `(x, y)` | MPAS-A |
| `dustin` | Dust tracers handed to other chemistry (used if `chem_opt` is 2 or 11) | ug kg-1 | `(x, y, 5)` | MPAS-A |
| `chem_opt` | WRF chemistry option. 2 and 11 write the update to `dustin`. Default 0. | 1 | scalar | constant |
| `ustar` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `znt` | Roughness length | m | `(x, y)` | MPAS-A |
| `snowh` | Snow depth | m | `(x, y)` | MPAS-A |
| `erod` | Erodibility (3 Ginoux classes, summed inside) | 1 | `(x, y, 3)` | external dataset |
| `clay` | Soil clay mass fraction | 1 | `(x, y)` | external dataset |
| `sand` | Soil sand mass fraction | 1 | `(x, y)` | external dataset |
| `oro` | Land-water mask (land if `oro < 1.5`) | 1 | `(x, y)` | MPAS-A |
| `isltyp` | Soil type index | 1 | `(x, y)` | external dataset |
| `dust_veg` | Vegetation option: 0 Ginoux mask, 1 vegetation fraction, 2 LAI mask | 1 | scalar | constant |
| `vegfra` | Vegetation fraction, needed if `dust_veg=1` | % | `(x, y)` | external dataset |
| `lai_vegmask` | MODIS vegetation mask, needed if `dust_veg=2` | 1 | `(x, y)` | external dataset |
| `dust_smois` | Moisture option (1 selects the RUC/PX dry limit) | 1 | scalar | constant |
| `sf_surface_physics` | WRF land surface model number (2 Noah, 3 RUC, 7 PX) | 1 | scalar | constant |
| `dust_dsr` | Option 1 selects the DRI erodibility | 1 | scalar | constant |
| `erod_dri` | DRI erodibility, needed if `dust_dsr=1` | 1 | `(x, y, 3)` | external dataset |
| `dust_soils` | Option 1 selects the NGA soil texture | 1 | scalar | constant |
| `clay_nga` | NGA clay fraction, needed if `dust_soils=1` | 1 | `(x, y)` | external dataset |
| `sand_nga` | NGA sand fraction, needed if `dust_soils=1` | 1 | `(x, y)` | external dataset |
| `settings` | Tunable source coefficients. See `tunable_names(AfwaSettings)`. | - | - | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the keys of the dict that `afwa_dust_driver` returns. The inputs are not changed.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `flux` | Dust emission flux per bin | kg m-2 s-1 | `(x, y, 5)` | MICM rate parameter |
| `dust` | Dust tracers after the emission (all levels kept, unchanged if `chem_opt` is 2 or 11) | ug kg-1 | `(x, level, y, 5)` | MPAS-A state |
| `dustin` | Updated lowest-level tracers if `chem_opt` is 2 or 11 | ug kg-1 | `(x, y, 5)` | MPAS-A state |
| `emis_dust` | Accumulated emitted dust per bin | kg m-2 | `(x, y, 5)` | diagnostic |
| `tot_edust` | Accumulated emitted dust, all bins | kg m-2 | `(x, y)` | diagnostic |
| `tot_dust` | Dust mass concentration, summed over bins | ug m-3 | `(x, level, y)` | diagnostic |
| `vis_dust` | Visibility | m | `(x, level, y)` | diagnostic |
| `afwa_dustloft` | Dust lofting potential (-99 over sea or where the erodibility is not positive) | m s-1 | `(x, y)` | diagnostic |

The tendency of the lowest layer is dC/dt = flux / (rho_air * dz).

## Usage

```python
import numpy as np
from quacs.dust.afwa import afwa_dust_driver

one = np.ones((1, 1))
level = np.ones((1, 1, 1))
out = afwa_dust_driver(
    dust=np.zeros((1, 1, 1, 5)),
    rho_phy=1.1 * level,
    dz8w=10.0 * level,
    smois=0.06 * level,
    dt=60.0,
    u10=8.0 * one,
    v10=0.0 * one,
    emis_dust=np.zeros((1, 1, 5)),
    tot_edust=np.zeros((1, 1)),
    dustin=np.zeros((1, 1, 5)),
    ustar=0.7 * one,
    znt=1e-3 * one,
    snowh=0.0 * one,
    erod=np.full((1, 1, 3), 0.2),
    clay=0.1 * one,
    sand=0.6 * one,
    oro=one,
    isltyp=3 * one,
)
F_p = out["flux"]  # kg m-2 s-1, shape (1, 1, 5)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the saltation thresholds and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 5 bin fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/afwa
```

## References

- LeGrand, S. L., et al. (2019). The AFWA dust emission scheme for the GOCART aerosol model in WRF-Chem v3.8.1. Geosci. Model Dev., 12, 131-166. https://doi.org/10.5194/gmd-12-131-2019
- Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme. J. Geophys. Res.-Atmos., 100, 16415-16430. https://doi.org/10.1029/95JD00690
- Iversen, J. D. and White, B. R. (1982). Saltation threshold on Earth, Mars and Venus. Sedimentology, 29, 111-119. https://doi.org/10.1111/j.1365-3091.1982.tb01713.x
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157. https://doi.org/10.1007/s00585-999-0149-7
- Kok, J. F. (2011). A scaling theory for the size distribution of emitted dust aerosols suggests climate models underestimate the size of the global dust cycle. PNAS, 108, 1016-1021. https://doi.org/10.1073/pnas.1014798108
