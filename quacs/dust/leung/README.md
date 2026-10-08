# Leung et al. (2023, 2024) dust emission

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Sue Park |
| **Original source** | CESM2/CLM5 `DustEmisLeung2023.F90` and `SoilStateInitTimeConstMod.F90`, [ESCOMP/CTSM](https://github.com/ESCOMP/CTSM) |

This scheme computes the vertical flux of mineral dust in 4 transport bins.
It follows the process-based and scale-aware scheme of Leung et al. (2023, 2024), as CTSM implements it.
It starts from the K14 flux equation (Kok et al., 2014b) and adds a drag partition for rock and vegetation and an intermittency factor (Comola et al., 2019).
The flux uses the impact threshold, so the scheme gives emission at lower winds in marginal source regions.
A time-invariant correction map `Kc_map` scales the flux, and lognormal source modes split it into the bins 0.1–1, 1–2.5, 2.5–5, and 5–10 um.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, radiative impacts |
| Scheme specifics | Leung et al. (2023, 2024) dust emission scheme, CESM2/CLM5 implementation. The dry threshold follows Shao & Lu (2000). The soil moisture threshold follows Fécan et al. (1999) with the Zender (2003a) form. The intermittency follows Comola et al. (2019). |
| External datasets | Clay fraction, saturated soil water content, aeolian roughness length `z0a` (Prigent et al., 2005), rock and vegetation area fractions, upscaling correction map `Kc_map` (Leung et al., 2024), land type (soil or crop mask). These are static fields that must be mapped to the MPAS grid. |
| Information from MPAS-A | Friction velocity, air density at the surface, volumetric soil water, liquid and frozen soil water, LAI, SAI, snow cover fraction, Obukhov length |
| Other inputs | Source-to-sink overlap matrix `ovr_src_snk_mss`, from `precompute_constants()`. Module constants: median soil diameter Dp = 127 um, VAI threshold 1.0, impact threshold ratio B_it = 0.82. |
| Information to MPAS-A | Dust emission flux `F_p` [x, y, 4 bins]. It must be added to the dust species of the host model. |
| Considerations | The 4 transport bins must be mapped to the dust species or bins of the host model. The intermittency factor needs the Obukhov length from the surface layer scheme. Dust has dry deposition. |
| Equations and references | u*s = u* F_eff (drag partition); eta = 1 - P_ft + alpha (P_ft - P_it); F_d = eta C_tune C_d f_bare f_clay' liqfrac rho (u*s^2 - u*it^2) / u*it (u*s / u*it)^kappa for u*s > u*it; F_p[n] = sum_m ovr_src_snk_mss[m, n] F_d Kc_map |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | sectional (4 bins) |
| Grid | 0-D |
| Interdependence | needs the Obukhov length from the surface layer scheme and the land surface fields from the host model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `dust_emission`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `fv` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `forc_rho` | Air density at the surface | kg m-3 | `(x, y)` | MPAS-A |
| `h2osoi_vol` | Volumetric soil water content, top layer | m3 m-3 | `(x, y)` | MPAS-A |
| `h2osoi_liq` | Liquid soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `h2osoi_ice` | Frozen soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `watsat` | Saturated volumetric soil water content | m3 m-3 | `(x, y)` | external dataset |
| `tlai` | One-sided leaf area index | m2 m-2 | `(x, y)` | MPAS-A |
| `tsai` | One-sided stem area index | m2 m-2 | `(x, y)` | MPAS-A |
| `frac_sno` | Snow cover fraction | 1 | `(x, y)` | MPAS-A |
| `obu` | Obukhov length | m | `(x, y)` | MPAS-A |
| `fclay` | Clay mass fraction | 1 | `(x, y)` | external dataset |
| `z0a` | Aeolian roughness length | m | `(x, y)` | external dataset |
| `frac_rock` | Rock area fraction A_r | 1 | `(x, y)` | external dataset |
| `frac_veg` | Vegetation area fraction A_v | 1 | `(x, y)` | external dataset |
| `Kc_map` | Upscaling correction factor | 1 | `(x, y)` | external dataset |
| `is_soil` | True for soil or crop land types | bool | `(x, y)` | external dataset |
| `ovr_src_snk_mss` | Source-to-sink mass overlap matrix (default from `precompute_constants()`) | 1 | `(3, 4)` | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `F_p` | Vertical dust emission flux in each transport bin (return value of `dust_emission`) | kg m-2 s-1 | `(x, y, 4)` | MPAS-A state |

The box model converts the flux to a tendency of the lowest layer: dC/dt = F_p / dz.

## Usage

```python
import numpy as np
from quacs.dust.leung import dust_emission

one = np.ones((1, 1))
F_p = dust_emission(
    fv=0.6 * one,
    forc_rho=1.2 * one,
    h2osoi_vol=0.05 * one,
    h2osoi_liq=5.0 * one,
    h2osoi_ice=0.0 * one,
    watsat=0.4 * one,
    tlai=0.1 * one,
    tsai=0.05 * one,
    frac_sno=0.0 * one,
    obu=-100.0 * one,
    fclay=0.1 * one,
    z0a=1.0e-4 * one,
    frac_rock=0.9 * one,
    frac_veg=0.1 * one,
    Kc_map=one,
    is_soil=np.ones((1, 1), dtype=bool),
)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the threshold, drag partition, and intermittency diagnostics and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 4 bin fluxes to a MUSICA MICM box model as `EMIS` rates |

## Tests

```bash
pytest quacs/dust/leung
```

## References

- Leung, D. M., et al. (2023). A new process-based and scale-aware desert dust emission scheme for global climate models – Part I: Description and evaluation against inverse modeling emissions. Atmos. Chem. Phys., 23, 6487–6523. https://doi.org/10.5194/acp-23-6487-2023
- Leung, D. M., et al. (2024). A new process-based and scale-aware desert dust emission scheme for global climate models – Part II: Evaluation in the Community Earth System Model version 2 (CESM2). Atmos. Chem. Phys., 24, 2287–2318. https://doi.org/10.5194/acp-24-2287-2024
- Kok, J. F., et al. (2014b). An improved dust emission model – Part 2: Evaluation in the Community Earth System Model, with implications for the use of dust source functions. Atmos. Chem. Phys., 14, 13043–13061. https://doi.org/10.5194/acp-14-13043-2014
- Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold friction velocity. J. Geophys. Res.-Atmos., 105, 22437–22443. https://doi.org/10.1029/2000JD900304
- Comola, F., Kok, J. F., Chamecki, M., and Martin, R. L. (2019). The intermittency of wind-driven sand transport. Geophys. Res. Lett., 46, 13430–13440. https://doi.org/10.1029/2019GL085739
- Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme. J. Geophys. Res., 100, 16415–16430. https://doi.org/10.1029/95JD00690
- Okin, G. S. (2008). A new model of wind erosion in the presence of vegetation. J. Geophys. Res.-Earth, 113. https://doi.org/10.1029/2007JF000758
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149–157. https://doi.org/10.1007/s00585-999-0149-7
