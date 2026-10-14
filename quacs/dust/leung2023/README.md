# Leung et al. (2023) dust emission (CESM2.1, CLM5 code)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Chayan Roychoudhury |
| **Original source** | `DUSTMod.F90` and `SoilStateInitTimeConstMod.F90` of the CESM2.1 fork `dmleung/dustmodelcesm2023`, [GitHub](https://github.com/dmleung/dustmodelcesm2023), commit `90771fefa1a36fb942eed6092a18e93896d01a4e`, archived as Zenodo record 10621844 version 1.0 (https://doi.org/10.5281/zenodo.10621844) |

This scheme computes the vertical flux of mineral dust in 4 size bins.
It follows the process-based and scale-aware scheme of Leung (2023), as the Zenodo code implements it.
The flux uses the Kok (2014) equation. It has an impact threshold, a drag partition for rock and vegetation, and an intermittency factor.
The driver adds the CTSM soil and patch steps and a coastal factor. It also adds the split into 4 bins and the dry deposition velocity.

The code follows the Zenodo Fortran, not the paper, where the two differ.

Checked against the compiled pinned Fortran: 12158 test cases, no differences at 8-byte precision.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms |
| Scheme specifics | Leung et al. (2023) dust emission, CESM2.1 code. The dry threshold follows Shao and Lu (2000). The moisture factor follows Fécan et al. (1999). The intermittency follows Comola et al. (2019), with atmospheric stability fixed to 0. The median soil diameter is 130 um in the code. |
| External datasets | Sand and clay percent, organic matter density, rock, vegetation and sparse-vegetation area fractions, rock drag factor, land fraction. These are static fields. Map them to the MPAS grid. |
| Information from MPAS-A | Friction velocity, air density, volumetric soil water, liquid and frozen soil water, LAI, SAI, snow cover fraction, surface pressure, air temperature, land unit and patch data |
| Other inputs | Aerodynamic resistance `ram1`. The setup dictionary from `init_dust_vars()` holds the bin overlap matrix, bin diameters and Stokes correction. |
| Information to MPAS-A | Dust emission flux in 4 bins. Add it to the dust species of the host model. The dry deposition velocity per bin is also given. |
| Considerations | Map the 4 bins (0.1 to 1, 1 to 2.5, 2.5 to 5, 5 to 10 um) to the dust species or bins of the host model. The land unit and patch inputs follow the CTSM data structure. The code uses `g = 9.80616`. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | sectional (4 bins) |
| Grid | 0-D |
| Interdependence | needs land surface fields (LAI, soil water, patches) from the host model and `ram1` from the surface layer scheme |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `leung_dust_driver`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `ustar` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `rhoa` | Air density at the surface | kg m-3 | `(x, y)` | MPAS-A |
| `h2osoi_vol` | Volumetric soil water, top layer | m3 m-3 | `(x, y)` | MPAS-A |
| `tlai` | One-sided leaf area index | m2 m-2 | `(x, y)` | MPAS-A |
| `tsai` | One-sided stem area index | m2 m-2 | `(x, y)` | MPAS-A |
| `frac_sno` | Snow cover fraction | 1 | `(x, y)` | MPAS-A |
| `pct_sand` | Sand percent | % | `(x, y)` | external dataset |
| `pct_clay` | Clay percent | % | `(x, y)` | external dataset |
| `organic` | Organic matter density | kg m-3 | `(x, y)` | external dataset |
| `h2osoi_liq` | Liquid soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `h2osoi_ice` | Frozen soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `land_fraction` | Land fraction of the cell (below 0.5 gets the coastal factor) | 1 | `(x, y)` | external dataset |
| `landunit_type` | Land unit type (1 soil, 2 crop emit) | 1 | `(x, y)` | MPAS-A |
| `patch_active` | True for active patches | bool | `(x, y)` | MPAS-A |
| `patch_wtlunit` | Patch weight on its land unit | 1 | `(x, y)` | MPAS-A |
| `rockfrc` | Rock area fraction A_r | 1 | `(x, y)` | external dataset |
| `vegefrc` | Shrub and crop area fraction A_v | 1 | `(x, y)` | external dataset |
| `sparfrc` | Sparse vegetation area fraction A_s | 1 | `(x, y)` | external dataset |
| `roughfct` | Rock drag partition factor F_eff | 1 | `(x, y)` | external dataset |
| `forc_pbot` | Surface pressure | Pa | `(x, y)` | MPAS-A |
| `forc_t` | Air temperature | K | `(x, y)` | MPAS-A |
| `ram1` | Aerodynamic resistance | s m-1 | `(x, y)` | other parameterization |
| `settings` | Tunable source coefficients. See `tunable_names(LeungSettings)`. | - | - | constant |
| `lake` | True for lake patches (flux is the fill value) | bool | `(x, y)` | MPAS-A |
| `landunit_id` | Land unit number of each patch, optional | 1 | `(x, y)` | MPAS-A |
| `pct_sand_layers` | Sand percent of all soil layers, optional input check | % | `(x, y, nlev)` | external dataset |
| `pct_clay_layers` | Clay percent of all soil layers, optional input check | % | `(x, y, nlev)` | external dataset |
| `setup` | Output of `init_dust_vars()` (made from the settings if not given) | - | dict | constant |
| `initial_value` | Fill value for lake patches and inactive patches | - | scalar | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.
These are the keys of the dictionary that `leung_dust_driver` returns.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `flx_mss_vrt_dst` | Vertical dust flux in each bin | kg m-2 s-1 | `(x, y, 4)` | MICM rate parameter |
| `flx_mss_vrt_dst_tot` | Sum of the 4 bins | kg m-2 s-1 | `(x, y)` | diagnostic |
| `flux_before_split` | Flux before the split into bins | kg m-2 s-1 | `(x, y)` | diagnostic |
| `watsat` | Saturated volumetric soil water | m3 m-3 | `(x, y)` | diagnostic |
| `gwc` | Gravimetric soil water | kg kg-1 | `(x, y)` | diagnostic |
| `liqfrac` | Liquid fraction of soil water | 1 | `(x, y)` | diagnostic |
| `tlai_lu` | Vegetation area index averaged over the land unit | m2 m-2 | `(x, y)` | diagnostic |
| `gwc_thr` | Threshold gravimetric soil water | kg kg-1 | `(x, y)` | diagnostic |
| `f_clay` | Clay fraction used in the flux (capped at 0.20) | 1 | `(x, y)` | diagnostic |
| `u_ft0` | Dry fluid threshold friction velocity | m s-1 | `(x, y)` | diagnostic |
| `u_ft` | Wet fluid threshold friction velocity | m s-1 | `(x, y)` | diagnostic |
| `u_it` | Impact threshold friction velocity | m s-1 | `(x, y)` | diagnostic |
| `u_st` | Standardized threshold friction velocity | m s-1 | `(x, y)` | diagnostic |
| `Cd` | Dust emission coefficient | 1 | `(x, y)` | diagnostic |
| `frag_expt` | Fragmentation exponent | 1 | `(x, y)` | diagnostic |
| `u_s` | Friction velocity after the drag partition | m s-1 | `(x, y)` | diagnostic |
| `lnd_frc_mbl` | Bare (mobile) land fraction | 1 | `(x, y)` | diagnostic |
| `fm` | Soil moisture factor on the threshold | 1 | `(x, y)` | diagnostic |
| `eta` | Intermittency factor | 1 | `(x, y)` | diagnostic |
| `f_eff` | Drag partition factor | 1 | `(x, y)` | diagnostic |
| `ovr_src_snk_mss` | Source-to-bin mass overlap matrix | 1 | `(3, 4)` | diagnostic |
| `dmt_vwr` | Volume-weighted diameter of each bin | m | `(4,)` | diagnostic |
| `stk_crc` | Stokes correction of each bin | 1 | `(4,)` | diagnostic |
| `dns_aer` | Dust particle density | kg m-3 | scalar | diagnostic |
| `vlc_trb` | Turbulent dry deposition velocity | m s-1 | `(x, y, 4)` | MPAS-A state |
| `vlc_grv` | Gravitational settling velocity | m s-1 | `(x, y, 4)` | diagnostic |
| `slp_crc` | Cunningham slip correction | 1 | `(x, y, 4)` | diagnostic |
| `rss_lmn` | Quasi-laminar surface resistance | s m-1 | `(x, y, 4)` | diagnostic |

The box model converts the flux to a tendency of the lowest layer: dC/dt = flux / dz.

## Usage

```python
import numpy as np
from quacs.dust.leung2023 import leung_dust_driver

one = np.ones((1, 1))
result = leung_dust_driver(
    ustar=0.7 * one,
    rhoa=1.1 * one,
    h2osoi_vol=0.02 * one,
    tlai=0.06 * one,
    tsai=0.04 * one,
    frac_sno=0.0 * one,
    pct_sand=65.0 * one,
    pct_clay=10.0 * one,
    organic=0.0 * one,
    h2osoi_liq=1.0 * one,
    h2osoi_ice=0.0 * one,
    land_fraction=one,
    landunit_type=np.array([[1]]),
    patch_active=np.array([[True]]),
    patch_wtlunit=one,
    rockfrc=0.5 * one,
    vegefrc=0.1 * one,
    sparfrc=0.1 * one,
    roughfct=0.7 * one,
    forc_pbot=1.0e5 * one,
    forc_t=295.0 * one,
    ram1=100.0 * one,
)
flux = result["flx_mss_vrt_dst"]
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the threshold, drag partition and intermittency values and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 4 bin fluxes to a MUSICA MICM box model as emission rates |

## Tests

```bash
pytest quacs/dust/leung2023
```

## References

- Leung, D. M., et al. (2023). A new process-based and scale-aware desert dust emission scheme for global climate models. Part I: Description and evaluation against inverse modeling emissions. ACP, 23, 6487-6523. https://doi.org/10.5194/acp-23-6487-2023
- Leung, D. M., et al. (2024). A new process-based and scale-aware desert dust emission scheme for global climate models. Part II: Evaluation in the Community Earth System Model version 2 (CESM2). ACP, 24, 2287-2318. https://doi.org/10.5194/acp-24-2287-2024
- Kok, J. F., et al. (2014). An improved dust emission model. Part 1: Model description and comparison against measurements. ACP, 14, 13023-13041. https://doi.org/10.5194/acp-14-13023-2014
- Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold friction velocity. J. Geophys. Res.-Atmos., 105, 22437-22443. https://doi.org/10.1029/2000JD900304
- Comola, F., et al. (2019). The intermittency of wind-driven sand transport. GRL, 46, 13430-13440. https://doi.org/10.1029/2019GL085739
- Okin, G. S. (2008). A new model of wind erosion in the presence of vegetation. J. Geophys. Res.-Earth, 113. https://doi.org/10.1029/2007JF000758
- Pierre, C., et al. (2014). Comparing drag partition schemes over a herbaceous Sahelian rangeland. J. Geophys. Res.-Earth, 119, 2291-2313. https://doi.org/10.1002/2014JF003177
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149-157. https://doi.org/10.1007/s00585-999-0149-7
