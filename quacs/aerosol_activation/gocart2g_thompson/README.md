# GOCART-2G aerosol activation for Thompson microphysics

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Aerosol activation and nucleation |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Aerosol-cloud interaction |
| Scheme specifics | GOCART-2G scheme (Collow et al., 2024); Thompson aerosol-aware scheme (Thompson et al., 2014) |
| External datasets | Aerosol initial and boundary conditions, emissions |
| Information from MPAS-A | T, P, W, Qv |
| Other inputs | Hygroscopicity (kappa) and mean number radius (r_m) |
| Information to MPAS-A | Water-friendly aerosol (WFA) and ice-friendly aerosol (IFA) number concentration [x, y, z] |
| Considerations | Derive r_m and kappa from GOCART-2G assumptions, not from fixed values. Other processes, such as coagulation, change the CCN size. |
| Equations and references | CCN activation lookup table from a parcel model [T, W, Nt, kappa, r_m] (Eidhammer et al., 2009) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | bulk |
| Grid | 0-D |
| Interdependence | needs GOCART-2G and Thompson microphysics |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `w` | Vertical velocity | m s-1 | to be defined | MPAS-A |
| `qv` | Water vapor mixing ratio | kg kg-1 | to be defined | MPAS-A |
| `aerosol_mass` | GOCART-2G aerosol mass per species | kg kg-1 | to be defined | MPAS-A |
| `kappa` | Hygroscopicity | 1 | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `nwfa` | Water-friendly aerosol number | kg-1 | to be defined | MPAS-A state |
| `nifa` | Ice-friendly aerosol number | kg-1 | to be defined | MPAS-A state |

## References

- Thompson, G. and Eidhammer, T. (2014). A study of aerosol impacts on clouds and precipitation development in a large winter cyclone. J. Atmos. Sci., 71, 3636-3658. https://doi.org/10.1175/JAS-D-13-0305.1
- Eidhammer, T., DeMott, P. J., and Kreidenweis, S. M. (2009). A comparison of heterogeneous ice nucleation parameterizations using a parcel model framework. J. Geophys. Res., 114, D06202. https://doi.org/10.1029/2008JD011095
- Collow, A., et al. (2024). GOCART-2G. (Full citation to be added.)
