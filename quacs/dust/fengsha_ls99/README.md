# FENGSHA dust emission, Lu and Shao (1999) form (CMAQ)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | CMAQ `CCTM/src/emis/emis/DUST_EMIS.F` (`get_dust_emis`, `dust_hflux`) and `lus_data_module.F`, [USEPA/CMAQ](https://github.com/USEPA/CMAQ), commit `9bd3734176479c2e49139fea98e1d5e8a16170e3` (no release tag) |

This scheme computes the vertical flux of mineral dust in 2 size modes (accumulation and coarse).
It follows the FENGSHA dust emission of CMAQ.
The horizontal saltation flux comes from the 10-m wind, the soil texture and the land use.
A drag partition for vegetation and a soil moisture correction reduce it.
Fixed sandblasting factors turn it into a vertical flux, and fixed mass fractions split it into the 2 modes.
It is a different scheme from the RRFS-SD FENGSHA in `quacs/dust/fengsha_rrfs`.

Checked against the compiled pinned Fortran: 1596 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, air quality |
| Scheme specifics | CMAQ FENGSHA dust emission. The threshold follows Shao and Lu (2000). The saltation flux follows White (1979). The drag partition uses land-use classes and a seasonal vegetation height. The soil table has 16 CMAQ soil types. Emission is zero over water, with rain above 0.0254 m, and with snow cover of 0.001 or more. |
| External datasets | Soil type index, vegetation fraction, land-use fractions for one of 5 land-use schemes (`usgs24`, `modis`, `modis_noah`, `nlcd_modis`, `nlcd40`). Optional soil texture (clay, coarse sand, fine-medium sand). These are static or monthly fields. Map them to the MPAS grid. |
| Information from MPAS-A | 10-m wind speed, volumetric soil moisture, air density at the surface, snow cover, rain, land-water mask, layer height, Jacobian field `jacobm`, day of year |
| Other inputs | Soil texture table, sandblasting factors, seasonal vegetation height, mode fractions and mode densities are in the settings object. |
| Information to MPAS-A | Dust emission rate in 2 modes, in g m-3 s-1. A factor `layer_height * jacobm * 1e-3` gives a flux in kg m-2 s-1. Add it to the dust species of the host model. |
| Considerations | Map the 2 modes to the dust species or bins of the host model. The code uses `g = 9.80622`. The code remaps the raw soil type index to the CMAQ table unless `soil_type_is_premapped` is true. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | modal (2 modes) |
| Grid | 0-D |
| Interdependence | needs land-use fractions and soil fields from the host model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `fengsha_dust_driver`, the CMAQ driver of the scheme. It takes the lowest air level of `dens`, and picks the desert classes from the full land-use fractions `luf`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `u10m` | 10-m wind speed | m s-1 | `(x, y)` | MPAS-A |
| `soim_vol` | Volumetric soil moisture | m3 m-3 | `(x, y)` | MPAS-A |
| `dens` | Air density (the code reads level 1) | kg m-3 | `(x, level, y)` | MPAS-A |
| `snocov` | Snow cover fraction (the code sets negative values to 0 and warns) | 1 | `(x, y)` | MPAS-A |
| `rn` | Non-convective rain | m | `(x, y)` | MPAS-A |
| `rc` | Convective rain | m | `(x, y)` | MPAS-A |
| `vegfrac` | Vegetation fraction | 1 | `(x, y)` | external dataset |
| `soil_type` | Soil type index (raw index, 1 to 16 after remap) | 1 | `(x, y)` | external dataset |
| `luf` | Full land-use fractions of the scheme `scheme` (24 columns for `usgs24`) | 1 | `(x, y, ncat)` | external dataset |
| `jdate` | Date as YYYYDDD (the day of year selects the seasonal vegetation height) | - | scalar | MPAS-A |
| `lwmask` | Land-water mask (above 0 is land) | 1 | `(x, y)` | MPAS-A |
| `layer_height` | Height of the lowest layer | m | `(x, y)` | MPAS-A |
| `jacobm` | Host Jacobian field. The driver forms `rjacm = 1 / jacobm` (ASX_DATA_MOD.F:878-882). With `mpas=True` the driver uses `rjacm = 1` and ignores `jacobm`, which may be None. | 1 | `(x, y)` | MPAS-A |
| `wmax_last` | Moisture limit saved from the previous call. Default None means zero (first step). Pass the returned `wmax` on later steps. | % | `(x, y)` | MPAS-A |
| `sd_ep_last` | Erodible potential saved from the previous call. Default None means zero (first step). Pass the returned `sd_ep` on later steps. | 1 | `(x, y)` | MPAS-A |
| `scheme` | Land-use scheme name (`usgs24`, `modis`, `modis_noah`, `nlcd_modis`, `nlcd40`). Default `usgs24`. | - | string | constant |
| `soil_type_is_premapped` | True if `soil_type` is already a CMAQ table row. Default False. | bool | scalar | constant |
| `clay` | Clay fraction, optional soil texture | 1 | `(x, y)` | external dataset |
| `csand` | Coarse sand fraction, optional soil texture | 1 | `(x, y)` | external dataset |
| `fmsand` | Fine-medium sand fraction, optional soil texture | 1 | `(x, y)` | external dataset |
| `diagnostic_totals` | Previous sum of the diagnostics, in rate-seconds. Optional. | - | `(x, y, n)` | MPAS-A |
| `sync_seconds` | Seconds added to the diagnostic sums per call. Default 1. | s | scalar | MPAS-A |
| `mpas` | True if `soim_vol` has a level axis (the code reads level 1). Default False. | bool | scalar | constant |
| `settings` | Tunable source coefficients. See `tunable_names(Ls99Settings)`. The species table `spcfac` (shape `(species, 2)`) defaults to the 28-row CMAQ `dust_spc` table, `AERO_DATA.F:505-532`. `DUST_SPC_NAMES` gives the row names. | - | - | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the main keys of the dict that `fengsha_dust_driver` returns. The dict has more diagnostic keys (`number_rates`, `surface_area_rates`, `qam`, `elus`, `ustr`, `kvh`, `fruf`, `soimt`, `fmoit`, `sd_ep`, `wmax`, `vegfrac`, `diagnostics`, `diagnostic_totals`).

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `mode_mass` | Dust emission rate per mode: accumulation, coarse. A rate per volume of the lowest layer. | g m-3 s-1 | `(x, y, 2)` | MICM rate parameter |
| `surface_rate` | Total dust emission flux of both modes | g m-2 s-1 | `(x, y)` | diagnostic |
| `dust_em` | Total dust emission rate of both modes | g m-3 s-1 | `(x, y)` | diagnostic |
| `species_rates` | `mode_mass` split into species with `spcfac` | g m-3 s-1 | `(x, y, 2, species)` | diagnostic |

For a flux per mode in kg m-2 s-1, multiply `mode_mass` by `layer_height` * `jacobm` (= `layer_height` / `rjacm`, and `rjacm` is 1 with `mpas=True`) and by 1e-3. Then the tendency of the lowest layer is dC/dt = F_p / (rho_air * dz).

## Usage

```python
import numpy as np
from quacs.dust.fengsha_ls99 import fengsha_dust_driver

one = np.ones((1, 1))
luf = np.zeros((1, 1, 24))
luf[0, 0, 18] = 0.8
out = fengsha_dust_driver(
    u10m=15.0 * one,
    soim_vol=0.02 * one,
    dens=1.2 * np.ones((1, 1, 1)),
    snocov=0.0 * one,
    rn=0.0 * one,
    rc=0.0 * one,
    vegfrac=0.05 * one,
    soil_type=np.array([[1]]),
    luf=luf,
    jdate=2023100,
    lwmask=one,
    layer_height=10.0 * one,
    jacobm=one,
    scheme="usgs24",
)
# g m-3 s-1 -> kg m-2 s-1: times layer_height (10 m) * jacobm (1) and 1e-3
F_p = out["mode_mass"] * 10.0 * 1.0 * 1e-3  # shape (1, 1, 2)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the flux for the land-use schemes and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 2 mode fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/fengsha_ls99
```

## References

The papers below come from the comments in the code. The module docstring cites none.

- Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold friction velocity. J. Geophys. Res.-Atmos., 105, 22437-22443. https://doi.org/10.1029/2000JD900304
- White, B. R. (1979). Soil transport by winds on Mars. J. Geophys. Res.-Sol. Ea., 84, 4643-4651. https://doi.org/10.1029/JB084iB09p04643
- Kok, J. F. (2011). A scaling theory for the size distribution of emitted dust aerosols suggests climate models underestimate the size of the global dust cycle. PNAS, 108, 1016-1021. https://doi.org/10.1073/pnas.1014798108
- Nabat, P., et al. (2012). Dust emission size distribution impact on aerosol budget and radiative forcing over the Mediterranean region. ACP, 12, 10545-10567. https://doi.org/10.5194/acp-12-10545-2012
- Darmenova, K., et al. (2009). Development of a physically based dust emission module within the Weather Research and Forecasting (WRF) model. J. Geophys. Res.-Atmos., 114. https://doi.org/10.1029/2008JD011236
- Xi, X., and Sokolik, I. N. (2015). Seasonal dynamics of threshold friction velocity and dust emission in Central Asia. J. Geophys. Res.-Atmos., 120, 1536-1564. https://doi.org/10.1002/2014JD022471.
  - Note: the code gives no title. This paper fits the three uses in the code (roughness density, monthly vegetation height and drag partition). The other 2015 paper by the same authors, 10.1002/2015JD024092, is about dust trends.
- Marticorena and Bergametti [JGR, 1997], saltation flux constant, as cited in `sources/CCTM__src__emis__emis__DUST_EMIS.F:1241`. The code gives no title.
- Chatenet, B., Marticorena, B., Gomes, L., and Bergametti, G. (1996). Assessing the microped size distributions of desert soils erodible by wind. Sedimentology, 43, 901-911. https://doi.org/10.1111/j.1365-3091.1996.tb01509.x
- Menut, L., et al. (2013). Impact of surface roughness and soil texture on mineral dust emission fluxes modeling. J. Geophys. Res.-Atmos., 118, 6505-6520. https://doi.org/10.1002/jgrd.50313
