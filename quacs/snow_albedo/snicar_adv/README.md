# SNICAR-ADv3 snow albedo

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Aerosol-Land interaction |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Aerosol-snow albedo interaction |
| Scheme specifics | Flanner et al. (2021), SNICAR-ADv3 solver |
| External datasets | Optics library for ice grains and each light-absorbing particle (LAP) species, 0.2-5 um. Refractive indices for ice and LAPs. Surface irradiance as a function of SZA. |
| Information from MPAS-A | SZA, spectral downward irradiance (direct and diffuse) |
| Other inputs | Snow properties from the land surface model: number of snow layers, thickness, snow depth, density, temperature of each layer. LAP deposition fluxes and mass mixing ratios. BC-snow mixing state. |
| Information to MPAS-A | Spectral snow albedo, absorbed solar flux in each snow layer |
| Considerations | Needs a land surface model. |
| Equations and references | Flanner et al. (2021) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | any |
| Grid | 1-D (column) |
| Interdependence | needs the land surface model |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `sza` | Solar zenith angle | degrees | to be defined | MPAS-A |
| `irradiance` | Spectral downward irradiance, direct and diffuse | W m-2 um-1 | to be defined | MPAS-A |
| `snow_layers` | Snow layer thickness, density, temperature | m, kg m-3, K | to be defined | other parameterization (LSM) |
| `lap_mass` | Light-absorbing particle mass in snow | kg kg-1 | to be defined | other parameterization |
| `optics_library` | Ice and LAP optical properties | various | to be defined | external dataset |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `snow_albedo` | Spectral snow albedo | 1 | to be defined | MPAS-A state |
| `absorbed_flux` | Absorbed solar flux in each snow layer | W m-2 | to be defined | MPAS-A state |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- Flanner et al., 2021, SNICAR-ADV 3 solver
