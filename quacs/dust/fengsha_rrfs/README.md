# FENGSHA dust emission, RRFS-SD form (Webb et al., 2020)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | RRFS-SD CCPP physics `physics/smoke_dust/dust_fengsha_mod.F90` and `dust_data_mod.F90`, [ufs-community/ccpp-physics](https://github.com/ufs-community/ccpp-physics), commit `242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7` (no release tag) |

This scheme computes the vertical flux of mineral dust in 5 size bins.
It is the NOAA/ARL FENGSHA scheme as the RRFS-SD (Rapid Refresh Forecast System, smoke and dust) uses it.
A horizontal saltation flux follows the Webb 2020 equation. It uses a drag partition and a soil moisture correction of the threshold.
A vertical-to-horizontal flux ratio turns it into a vertical flux. A lognormal distribution (Kok, 2011) splits it into 5 bins.

`rrfs_dust_step` also ports the RRFS host tracer steps. It takes the two RRFS dust tracers (dust and coarse PM). It builds the 5 bins and returns the two updated tracers.
It is not the CMAQ FENGSHA in `quacs/dust/fengsha_ls99`.

Checked against the compiled pinned Fortran: 576 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, air quality forecasting |
| Scheme specifics | FENGSHA of RRFS-SD. The flux is `(alpha / g) * ssm^gamma * rho_a * kvh * q`, with `q = max(0, R*u* - u_t*h) * (R*u* + u_t*h)^2`. The ratio `kvh` follows Marticorena and Bergametti (1995). The factor `h` is the Fécan et al. (1999) correction (default) or the Shao correction. The code converts soil water from volumetric to gravimetric. |
| External datasets | Clay, sand, sediment supply map `ssm`, drag partition `rdrag`, dry threshold friction velocity `uthrs`, soil type. These are static or monthly fields. Map them to the MPAS grid. The monthly driver `fengsha_rrfs_monthly_driver` selects the month. |
| Information from MPAS-A | Friction velocity, volumetric soil moisture, pressure and temperature of the lowest layer, layer pressure edges, cell area, soil temperature, roughness length, snow height, land type, time step, dust tracers |
| Other inputs | Scaling `alpha`, exponent `gamma`, `kvhmax`, moisture option, bin radii and tracer remap factors are in the settings object. |
| Information to MPAS-A | Dust emission flux in 5 bins and the updated dust tracers. Add them to the dust species of the host model. |
| Considerations | Map the 5 bins to the dust species or bins of the host model. The code uses `g = 9.80665`. Cells with frozen soil, snow, large roughness, non-erodible soil type or `rdrag <= 0` emit nothing. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | sectional (5 bins) |
| Grid | 0-D |
| Interdependence | needs friction velocity from the surface layer scheme and soil fields from the land surface model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `fengsha_rrfs_driver`, the port of the Fortran dust driver, which takes 5 bins.
`rrfs_dust_step` is the wrapper-level entry. It takes the host arrays (see its docstring). It takes two tracers, `dust` and `coarse_pm` (ug kg-1, shape `(x, y)`), in place of the 5 bins.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `ustar` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `slc` | Volumetric soil moisture, layer 2 | m3 m-3 | `(x, y)` | MPAS-A |
| `clay` | Clay mass fraction | 1 | `(x, y)` | external dataset |
| `sand` | Sand mass fraction | 1 | `(x, y)` | external dataset |
| `ssm` | Sediment supply map | 1 | `(x, y)` | external dataset |
| `rdrag` | Drag partition R (cells with `rdrag <= 0` emit nothing) | 1 | `(x, y)` | external dataset |
| `uthrs` | Dry threshold friction velocity | m s-1 | `(x, y)` | external dataset |
| `isltyp` | Integer soil type | 1 | `(x, y)` | external dataset |
| `stemp` | Soil temperature, layer 1 | K | `(x, y)` | MPAS-A |
| `znt` | Roughness length | m | `(x, y)` | MPAS-A |
| `snowh` | Snow height | m | `(x, y)` | MPAS-A |
| `xland` | Land type (land if below 1.5) | 1 | `(x, y)` | MPAS-A |
| `p` | Pressure of layer 1 | Pa | `(x, y)` | MPAS-A |
| `temperature` | Temperature of layer 1 | K | `(x, y)` | MPAS-A |
| `p8w1` | Pressure at the bottom edge of layer 1 | Pa | `(x, y)` | MPAS-A |
| `p8w2` | Pressure at the top edge of layer 1 | Pa | `(x, y)` | MPAS-A |
| `area` | Cell area | m2 | `(x, y)` | MPAS-A |
| `dt` | Time step | s | scalar | MPAS-A |
| `dust` | Initial dust tracers | ug kg-1 | `(x, y, 5)` | MPAS-A |
| `settings` | Tunable source coefficients. See `tunable_names(RrfsSettings)`. | - | - | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the fields of `FengshaRRFSResult`, the return value of `fengsha_rrfs_driver`.
`rrfs_dust_step` returns `RrfsStepResult` instead: `rate`, `bems`, `emis_dust`, and the two tracers `dust` and `coarse_pm` (`(x, y)`, ug kg-1). It floors them at `epsilc` (1e-12) and caps them at 200 and 5000 (W:589-590).

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `rate` | Dust emission flux per bin | kg m-2 s-1 | `(x, y, 5)` | MICM rate parameter |
| `bems` | Dust emission flux per bin (`rate` times 1e9) | ug m-2 s-1 | `(x, y, 5)` | diagnostic |
| `dust` | Updated dust tracers after the emission and the RRFS remap | ug kg-1 | `(x, y, 5)` | MPAS-A state |
| `emis_dust` | Sum of the diagnostic bins 1 to 4 of `bems` | ug m-2 s-1 | `(x, y)` | diagnostic |

The tendency of the lowest layer is dC/dt = rate / (rho_air * dz).

## Usage

```python
import numpy as np
from quacs.dust.fengsha_rrfs import rrfs_dust_step

one = np.ones((1, 1))
result = rrfs_dust_step(
    dust12m=np.tile([0.1, 1.0, 0.6, 1.0, 0.2], (1, 1, 12, 1)),  # clay, rdrag, sand, ssm, uthrs
    month=1,
    zorl=one,  # cm
    snow=0.0 * one,  # mm
    land=np.array([[1]]),
    smc=np.array([[[0.3, 0.05]]]),  # soil layers, layer 2 is used
    tslb=np.array([[[290.0]]]),  # soil layers, layer 1 is used
    pr3d=np.array([[[1.0e5, 0.99e5]]]),  # layer edges
    prl3d=np.array([[[1.0e5]]]),  # layer middles
    tk3d=np.array([[[300.0]]]),
    ustar=0.7 * one,
    garea=one,
    soiltyp=np.array([[2]]),
    dt=1.0,
    dust=0.0 * one,  # RRFS dust tracer, ug/kg
    coarse_pm=0.0 * one,  # RRFS coarse PM tracer, ug/kg
)
result.dust, result.coarse_pm  # updated tracers
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the threshold and moisture factor and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 5 bin fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/fengsha_rrfs
```

## References

- Webb et al. (2020), Eq. 9, as cited in `sources/physics__smoke_dust__dust_fengsha_mod.F90:483`. The code gives no title.
- Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme. J. Geophys. Res.-Atmos., 100, 16415-16430. https://doi.org/10.1029/95JD00690
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157. https://doi.org/10.1007/s00585-999-0149-7
- Kok, J. F. (2011). A scaling theory for the size distribution of emitted dust aerosols suggests climate models underestimate the size of the global dust cycle. PNAS, 108, 1016-1021. https://doi.org/10.1073/pnas.1014798108
- MacKinnon, D. J., et al. (2004). Comparison of aerodynamically and model-derived roughness lengths (zo) over diverse surfaces, central Mojave Desert, California, USA. Geomorphology, 63, 103-113. https://doi.org/10.1016/j.geomorph.2004.03.009
