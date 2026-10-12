# GOCART dust emission (Ginoux et al., 2001)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Dust Emission Scheme |
| **Contributors** | [chayanroyc](https://github.com/chayanroyc) |
| **Original source** | [MPAS-GOCART2G](https://github.com/PACE-DAAQ/MPAS-GOCART2G): `DU2G_GridCompMod.F90` (`DustEmissionGOCART2G_revised`) and `DU2G_instance.F90` |

This scheme computes the vertical flux of mineral dust from the surface in 5 size bins.
The flux increases with the cube of the 10-m wind speed above a threshold.
Soil moisture increases the threshold. The scheme sets the flux to zero over water, over lakes, and on wet soil.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, radiative impacts |
| Scheme specifics | GOCART dust emission scheme (Ginoux et al., 2001), as in GOCART-2G. Dry-soil threshold from Marticorena & Bergametti (1995). |
| External datasets | Erodibility `S`, land-water mask `oro`, lake fraction `frlake`, and the scale factor for each bin `s_p`. These are static fields that must be mapped to the MPAS grid. |
| Information from MPAS-A | 10-m wind speed `w10m`, soil moisture `gwet`, near-surface air density `rhoa` |
| Other inputs | Empirical constant `C` (set in the script or in a namelist). Particle radius and density for each bin. |
| Information to MPAS-A | Dust emission flux `F_p` [x, y, dust bin]. The host model adds it to the dust concentrations. |
| Considerations | The dust species must be in MPAS-A and mapped to its dust bins, or MICM must add them. Dust has dry deposition. Heterogeneous reactions on dust change ozone. |
| Equations and references | `u_thresh = max(0, u_t (1.2 + 0.2 log10(max(gwet, 1e-3))))`; `F_p = C S s_p (1 - frlake) w10m^2 max(w10m - u_thresh, 0)` (Ginoux et al., 2001, Eq. 1) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver (emission rate `F_p / dz`) |
| Aerosol representation | sectional (5 dust bins) |
| Grid | 0-D (surface cell) |
| Interdependence | standalone |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.

Inputs of `dust_emission(w10m, u_t, gwet, oro, frlake, S, s_p, C)`:

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `w10m` | 10-m wind speed | m s-1 | `(x, y)` | MPAS-A |
| `u_t` | Dry-soil threshold wind speed, from `threshold_velocity` | m s-1 | `(x, y)` or scalar | other parameterization (`threshold_velocity`) |
| `gwet` | Gravimetric soil moisture, in [0, 1] | 1 | `(x, y)` | MPAS-A |
| `oro` | Land-water mask (1.0 = land) | 1 | `(x, y)` | external dataset |
| `frlake` | Lake fraction, in [0, 1] | 1 | `(x, y)` | external dataset |
| `S` | Source erodibility, in [0, 1] | 1 | `(x, y)` | external dataset |
| `s_p` | Scale factor for each bin (default `SP_DEFAULT`) | 1 | `(n_bins,)` | external dataset |
| `C` | Empirical scale constant (default `C_DEFAULT = 0.088`) | kg s2 m-5 | scalar | constant |

Inputs of `threshold_velocity(radius, rhop, rhoa)`:

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `radius` | Effective particle radius for each bin (default `RADIUS_DEFAULT`) | m | `(n_bins,)` | constant |
| `rhop` | Particle density for each bin (default `RHOP_DEFAULT`) | kg m-3 | `(n_bins,)` | constant |
| `rhoa` | Near-surface air density | kg m-3 | `(x, y)` | MPAS-A |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `F_p` | Dust emission flux for each bin, from `dust_emission` | kg m-2 s-1 | `(x, y, n_bins)` | MPAS-A state (tendency `F_p / (rho_air dz)`) |
| `u_t` | Dry-soil threshold wind speed, from `threshold_velocity` | m s-1 | `(x, y, n_bins)` | diagnostic |
| `EMIS.dust_emis_bin_i` | Emission rate for each bin, from `emission_rates` (`F_p / dz`) | kg m-3 s-1 | `(1,)` for each bin | MICM rate parameter |

## Usage

```python
import numpy as np
from quacs.dust.ginoux import LAND, SP_DEFAULT, dust_emission

F_p = dust_emission(
    w10m=np.array([[12.0]]),
    u_t=np.array([[6.0]]),
    gwet=np.array([[0.1]]),
    oro=np.array([[LAND]]),
    frlake=np.array([[0.0]]),
    S=np.array([[0.5]]),
    s_p=SP_DEFAULT,
)
```

| Example | What it does |
|---|---|
| `examples/musica_box_model.py` | Runs the scheme in a MUSICA box model and compares the result with the closed-form solution |
| `examples/ginoux_box_model.ipynb` | Runs the box model and plots the dust mass for each bin and the flux as a function of wind speed |

## Tests

```bash
pytest quacs/dust/ginoux
```

## References

- Ginoux, P., et al. (2001). Sources and distributions of dust aerosols simulated with the GOCART model. J. Geophys. Res., 106, 20255-20273. https://doi.org/10.1029/2000JD000053
- Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme. J. Geophys. Res., 100, 16415-16430. https://doi.org/10.1029/95JD00690
