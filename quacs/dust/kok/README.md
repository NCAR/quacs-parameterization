# Kok et al. (2014) dust emission (K14)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Sue Park |
| **Original source** | CESM2/CLM5 `DUSTMod.F90` (K14 branch), [ESCOMP/CTSM](https://github.com/ESCOMP/CTSM) |

This scheme computes the total vertical flux of mineral dust from the surface.
It follows the K14 scheme of Kok et al. (2014a, 2014b), as CLM5 implements it.
The soil erodibility and the fragmentation exponent depend on the standardized threshold friction velocity, so the scheme needs no external dust source function.
The output is one total flux. The scheme does not split the flux into size bins.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, radiative impacts |
| Scheme specifics | K14 dust emission scheme (Kok et al., 2014a, 2014b), CLM5 implementation. The dry threshold follows Shao & Lu (2000). The soil moisture correction follows Fécan et al. (1999). |
| External datasets | Clay fraction and land-water mask (static fields that must be mapped to the MPAS grid) |
| Information from MPAS-A | Friction velocity, gravimetric soil moisture, air density at the surface, vegetation area index (LAI + SAI), snow depth |
| Other inputs | Tuning constant `C_tune` (default 0.05), median soil particle diameter `Dp` (default 75 um), particle density `rhop`, Fécan coefficient `a`, vegetation threshold `vai_thr` (default 0.3) |
| Information to MPAS-A | Total dust emission flux `F_d` [x, y]. It must be split into the dust species or modes of the host model. |
| Considerations | The scheme gives one total flux. The host model must split it into size bins or modes, for example with the brittle-fragmentation size distribution (Kok, 2011). The clay fraction in the flux is capped at 0.2. Dust has dry deposition. |
| Equations and references | u_st = u_ft sqrt(rho_a / rho_a0); C_d = C_d0 exp(-C_e (u_st - u_st0) / u_st0); kappa = C_kappa (u_st - u_st0) / u_st0; F_d = C_tune C_d f_bare f_clay rho_a (u*^2 - u_t^2) / u_st (u* / u_t)^kappa for u* > u_t |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | discrete state update |
| Aerosol representation | bulk (one total flux) |
| Grid | 0-D |
| Interdependence | standalone (needs the land surface fields from the host model) |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
These are the arguments of `dust_emission`.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `ustar` | Surface friction velocity | m s-1 | `(x, y)` | MPAS-A |
| `w` | Gravimetric soil moisture | kg kg-1 | `(x, y)` | MPAS-A |
| `rho_a` | Air density at the surface | kg m-3 | `(x, y)` | MPAS-A |
| `vai` | Vegetation area index (LAI + SAI) | m2 m-2 | `(x, y)` | MPAS-A |
| `snowdp` | Snow depth | m | `(x, y)` | MPAS-A |
| `f_clay` | Clay mass fraction | 1 | `(x, y)` | external dataset |
| `oro` | Land-water mask (1 = land) | 1 | `(x, y)` | external dataset |
| `Dp` | Median soil particle diameter (default 75e-6) | m | scalar | constant |
| `rhop` | Soil particle density (default 2650) | kg m-3 | scalar | constant |
| `a` | Fécan et al. (1999) tuning coefficient (default 1.0) | 1 | scalar | constant |
| `vai_thr` | VAI above which no dust is emitted (default 0.3) | m2 m-2 | scalar | constant |
| `C_tune` | Global tuning constant (default 0.05) | 1 | scalar | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `F_d` | Total vertical dust emission flux (return value of `dust_emission`) | kg m-2 s-1 | `(x, y)` | MPAS-A state |

The box model converts the flux to a tendency of the lowest layer: dx/dt = F_d / (rho_air dz).

## Usage

```python
import numpy as np
from quacs.dust.kok import dust_emission

F_d = dust_emission(
    ustar=np.array([[0.5]]),
    w=np.array([[0.02]]),
    rho_a=np.array([[1.2]]),
    vai=np.array([[0.05]]),
    snowdp=np.array([[0.0]]),
    f_clay=np.array([[0.15]]),
    oro=np.array([[1.0]]),
)
```

| Example | What it does |
|---|---|
| `examples/numpy_box_model.py` | Prints the threshold diagnostics and steps one cell forward with NumPy |
| `examples/musica_box_model.py` | Adds the flux to a MUSICA MICM box model as an `EMIS` rate |

## Tests

```bash
pytest quacs/dust/kok
```

## References

- Kok, J. F., et al. (2014a). An improved dust emission model – Part 1: Model description and comparison against measurements. Atmos. Chem. Phys., 14, 13023–13041. https://doi.org/10.5194/acp-14-13023-2014
- Kok, J. F., et al. (2014b). An improved dust emission model – Part 2: Evaluation in the Community Earth System Model, with implications for the use of dust source functions. Atmos. Chem. Phys., 14, 13043–13061. https://doi.org/10.5194/acp-14-13043-2014
- Shao, Y. and Lu, H. (2000). A simple expression for wind erosion threshold friction velocity. J. Geophys. Res., 105, 22437–22443. https://doi.org/10.1029/2000JD900304
- Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization of the increase of the aeolian erosion threshold wind friction velocity due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17, 149–157. https://doi.org/10.1007/s00585-999-0149-7
