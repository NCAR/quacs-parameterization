# Resolved-scale trace gas wet deposition (Neu et al., 2012)

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Resolved scale Trace Gas Wet Deposition |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Chemistry related to ozone, air quality |
| Scheme specifics | Neu et al. (2012). Code in WRF-Chem for MOZART chemistry schemes. |
| External datasets | None |
| Information from MPAS-A | p, T, air density, water vapor, cloud water, cloud ice, rain, snow, graupel, cloud fraction, precipitation production (rainprod), evaporation (evapprod), chemical species, time step, grid cell area, dz |
| Other inputs | Gas-aqueous fraction from Henry's law coefficients, retention fraction in freezing drops, molecular weight, effective Henry's law coefficient |
| Information to MPAS-A | Updated chemical species; HNO3 column change from scavenging (diagnostic) |
| Considerations | The scheme uses different equations for the all-liquid, mixed-phase, and all-ice regions of a storm. |
| Equations and references | Neu and Prather (2012) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | discrete state update |
| Aerosol representation | none |
| Grid | 1-D (column) |
| Interdependence | needs host microphysics |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `hydrometeors` | Cloud water, cloud ice, rain, snow, graupel | kg kg-1 | to be defined | MPAS-A |
| `rainprod` | Precipitation production | kg kg-1 s-1 | to be defined | MPAS-A |
| `evapprod` | Precipitation evaporation | kg kg-1 s-1 | to be defined | MPAS-A |
| `henry` | Henry's law coefficients | M atm-1 | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `species` | Updated gas-phase species | mol mol-1 | to be defined | MPAS-A state |
| `hno3_column_wetdep` | HNO3 column change from scavenging | kg m-2 s-1 | to be defined | diagnostic |

## References

- Neu, J. L. and Prather, M. J. (2012). Toward a more physical representation of precipitation scavenging in global chemistry models: cloud overlap and ice physics and their impact on tropospheric ozone. Atmos. Chem. Phys., 12, 3289-3310. https://doi.org/10.5194/acp-12-3289-2012
