# Aerosol optical properties (Ghan & Zaveri, 2007)

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Aerosol physics interaction |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Aerosol direct radiative effects, aerosol-radiation interactions |
| Scheme specifics | Ghan & Zaveri (2007) |
| External datasets | Mie scattering lookup tables for extinction and scattering efficiency and asymmetry, as a function of wet radius and refractive index |
| Information from MPAS-A | Temperature, pressure, humidity |
| Other inputs | Aerosol concentration for each species and mode, hygroscopicity (kappa) |
| Information to MPAS-A | AOD, SSA, and g for each radiation band |
| Considerations | Internal mixing in each mode. The scheme must compute the insoluble species, the refractive index, and the wet radius for each mode. The output goes to RRTMG in MPAS-A. |
| Equations and references | Ghan & Zaveri (2007) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | modal |
| Grid | 1-D (column) |
| Interdependence | needs the radiation scheme (RRTMG) |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `qv` | Specific humidity | kg kg-1 | to be defined | MPAS-A |
| `aerosol_mass` | Aerosol mass for each species and mode | kg kg-1 | to be defined | MPAS-A |
| `kappa` | Hygroscopicity for each species | 1 | to be defined | constant |
| `mie_tables` | Mie lookup tables | various | to be defined | external dataset |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `aod` | Aerosol optical depth for each band | 1 | to be defined | MPAS-A state (radiation) |
| `ssa` | Single-scattering albedo for each band | 1 | to be defined | MPAS-A state (radiation) |
| `g` | Asymmetry parameter for each band | 1 | to be defined | MPAS-A state (radiation) |

## References

- Ghan, S. J. and Zaveri, R. A. (2007). Parameterization of optical properties for hydrated internally mixed aerosol. J. Geophys. Res., 112, D10201. https://doi.org/10.1029/2006JD007927
