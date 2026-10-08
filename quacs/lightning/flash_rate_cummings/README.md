# Lightning flash rate for resolved convection

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Lightning flash rate using resolved convection |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Lightning-NOx effects on ozone, aerosols, aerosol-cloud interactions |
| Scheme specifics | Cummings et al. (2024). Information from WRF-Chem. |
| External datasets | None, unless a dataset gives the IC to CG ratio |
| Information from MPAS-A | Land or ocean cover, radar reflectivity, w, temperature, height, terrain height. Other parameterizations use wmax, updraft volume, qice, qsnow, qgraupel, pressure. |
| Other inputs | None |
| Information to MPAS-A | Intracloud flash rate IC_FR and cloud-to-ground flash rate CG_FR (flashes s-1) |
| Considerations | zkm is a bulk storm property, so the scheme must find a storm region. |
| Equations and references | total_FR = 3.44e-5 zkm^4.9 / 60 (land); 6.57e-6 zkm^4.9 / 60 (ocean); distributed to cells with refl > 20 dBZ; zkm = height above ground of the highest 20 dBZ level with T < 0 C |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | discrete state update |
| Aerosol representation | none |
| Grid | 1-D (column) |
| Interdependence | needs host microphysics (reflectivity) |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `land_mask` | Land or ocean flag | 1 | to be defined | MPAS-A |
| `reflectivity` | Radar reflectivity | dBZ | to be defined | MPAS-A |
| `w` | Vertical velocity | m s-1 | to be defined | MPAS-A |
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `height` | Height of each level | m | to be defined | MPAS-A |
| `terrain_height` | Terrain height | m | to be defined | MPAS-A |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `ic_flash_rate` | Intracloud flash rate | s-1 | to be defined | other parameterization (lightning NOx) |
| `cg_flash_rate` | Cloud-to-ground flash rate | s-1 | to be defined | other parameterization (lightning NOx) |

## References

- Cummings, K. A., et al. (2024). (Full citation to be added.)
