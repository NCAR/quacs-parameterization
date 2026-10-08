# DEAD dust emission (Zender et al., 2003)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Sue Park |
| **Original source** | CLM/CTSM `src/biogeochem/DustEmisZender2003.F90` and `src/biogeophys/SoilStateInitTimeConstMod.F90`, [ESCOMP/CTSM](https://github.com/ESCOMP/CTSM) |

This scheme computes the vertical flux of mineral dust in 4 transport bins.
It follows the Mineral Dust Entrainment and Deposition (DEAD) model of Zender et al. (2003a), as CLM4 implements it.
A horizontal saltation flux (White, 1979) with the Owen effect becomes a vertical flux through the sandblasting efficiency of Marticorena & Bergametti (1995).
Three lognormal source modes then split the flux into the bins 0.1–1, 1–2.5, 2.5–5, and 5–10 um.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, radiative impacts |
| Scheme specifics | DEAD dust emission scheme (Zender et al., 2003a), CLM4 implementation. The soil moisture correction follows Fécan et al. (1999) with the Zender tuning. |
| External datasets | Clay fraction, saturated soil water content, land type (soil or crop mask), basin erodibility factor (1 in CLM4). These are static fields that must be mapped to the MPAS grid. |
| Information from MPAS-A | Friction velocity, 10-m wind speed, air density at the surface, volumetric soil water, liquid and frozen soil water, LAI, SAI, snow cover fraction |
| Other inputs | Threshold friction velocity factor `tmp1` and the source-to-sink overlap matrix `ovr_src_snk_mss`. `precompute_constants()` computes both. |
| Information to MPAS-A | Dust emission flux `F_p` [x, y, 4 bins]. It must be added to the dust species of the host model. |
| Considerations | The 4 transport bins must be mapped to the dust species or bins of the host model. The clay fraction in the sandblasting efficiency is capped at 0.2. Dust has dry deposition. |
| Equations and references | u*_t = tmp1 / sqrt(rho) f_w; u*_s = u* + 0.003 (u10 - u10_t)^2; Q = c rho u*_s^3 / g (1 - u*_t/u*_s)(1 + u*_t/u*_s)^2 f_bare; F_tot = Q 100 10^(13.4 min(f_clay, 0.2) - 6); F_p[n] = sum_m ovr_src_snk_mss[m, n] F_tot |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | discrete state update |
| Aerosol representation | sectional (4 bins) |
| Grid | 0-D |
| Interdependence | standalone (needs the land surface fields from the host model) |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `dust_emission`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `fv` | Friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `u10` | 10-m wind speed | m s-1 | `(x, y)` | MPAS-A |
| `forc_rho` | Air density at the surface | kg m-3 | `(x, y)` | MPAS-A |
| `h2osoi_vol` | Volumetric soil water content, top layer | m3 m-3 | `(x, y)` | MPAS-A |
| `h2osoi_liq` | Liquid soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `h2osoi_ice` | Frozen soil water, top layer | kg m-2 | `(x, y)` | MPAS-A |
| `watsat` | Saturated volumetric soil water content | m3 m-3 | `(x, y)` | external dataset |
| `tlai` | One-sided leaf area index | m2 m-2 | `(x, y)` | MPAS-A |
| `tsai` | One-sided stem area index | m2 m-2 | `(x, y)` | MPAS-A |
| `frac_sno` | Snow cover fraction | 1 | `(x, y)` | MPAS-A |
| `fclay` | Clay mass fraction | 1 | `(x, y)` | external dataset |
| `mbl_bsn_fct` | Basin erodibility factor (1 in CLM4) | 1 | `(x, y)` | external dataset |
| `is_soil` | True for soil or crop land types | bool | `(x, y)` | external dataset |
| `tmp1` | Threshold friction velocity factor (default from `precompute_constants()`) | kg0.5 m-0.5 s-1 | scalar | constant |
| `ovr_src_snk_mss` | Source-to-sink mass overlap matrix (default from `precompute_constants()`) | 1 | `(3, 4)` | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `F_p` | Vertical dust emission flux in each transport bin (return value of `dust_emission`) | kg m-2 s-1 | `(x, y, 4)` | MPAS-A state |

The box model converts the flux to a tendency of the lowest layer: dC/dt = F_p / dz.

## Usage

```python
import numpy as np
from quacs.dust.zender import dust_emission

one = np.ones((1, 1))
F_p = dust_emission(
    fv=0.5 * one,
    u10=8.0 * one,
    forc_rho=1.2 * one,
    h2osoi_vol=0.05 * one,
    h2osoi_liq=5.0 * one,
    h2osoi_ice=0.0 * one,
    watsat=0.4 * one,
    tlai=0.0 * one,
    tsai=0.0 * one,
    frac_sno=0.0 * one,
    fclay=0.1 * one,
    mbl_bsn_fct=one,
    is_soil=np.ones((1, 1), dtype=bool),
)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the precomputed constants and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the 4 bin fluxes to a MUSICA MICM box model as `EMIS` rates. The forcing matches the Ginoux box model. |

## Tests

```bash
pytest quacs/dust/zender
```

## References

- Zender, C. S., Bian, H. S., and Newman, D. (2003a). Mineral Dust Entrainment and Deposition (DEAD) model: Description and 1990s dust climatology. J. Geophys. Res.-Atmos., 108, 4416. https://doi.org/10.1029/2002JD002775
- Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme. J. Geophys. Res., 100, 16415–16430. https://doi.org/10.1029/95JD00690
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149–157. https://doi.org/10.1007/s00585-999-0149-7
- Kok, J. F., et al. (2014). An improved dust emission model – Part 2: Evaluation in the Community Earth System Model, with implications for the use of dust source functions. Atmos. Chem. Phys., 14, 13043–13061. https://doi.org/10.5194/acp-14-13043-2014
