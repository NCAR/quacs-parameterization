# Resolved-scale aerosol wet deposition (MOSAIC)

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Resolved scale Aerosol Wet Deposition |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Aerosol life cycle, air quality, climate |
| Scheme specifics | PNNL code in WRF-Chem for the MOSAIC sectional scheme. It includes in-cloud scavenging and impaction (below-cloud) scavenging. |
| External datasets | None |
| Information from MPAS-A | p, T, air density, fractional cloud water, precipitation rate of rain, snow, and graupel, cloud fraction, aerosol mass and number concentrations |
| Other inputs | Aerosol density, bin sections, below-cloud scavenging coefficients |
| Information to MPAS-A | Updated aerosol concentrations; column change from scavenging (diagnostic) |
| Considerations | The scheme needs the vertical structure of the host model. Impaction scavenging starts at the model top and goes down to get the precipitation flux in each layer. |
| Equations and references | WRF-Chem module_mosaic_wetscav.F (to be added) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | sectional |
| Grid | 1-D (column) |
| Interdependence | needs host microphysics |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `air_density` | Air density | kg m-3 | to be defined | MPAS-A |
| `cloud_fraction` | Cloud fraction | 1 | to be defined | MPAS-A |
| `precip_rate` | Precipitation rate of rain, snow, graupel | kg m-2 s-1 | to be defined | MPAS-A |
| `aerosol` | Aerosol mass and number per bin | kg kg-1, kg-1 | to be defined | MPAS-A |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `aerosol` | Updated aerosol mass and number per bin | kg kg-1, kg-1 | to be defined | MPAS-A state |
| `column_wetdep` | Column change from scavenging | kg m-2 s-1 | to be defined | diagnostic |

## References

- Easter, R. C., et al. (2004). MIRAGE: Model description and evaluation of aerosols and trace gases. J. Geophys. Res., 109, D20210. https://doi.org/10.1029/2004JD004571
