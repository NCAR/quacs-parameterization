# Lightning flash rate for parameterized convection

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Lightning flash rate using parameterized convection |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Lightning-NOx effects on ozone, aerosols, aerosol-cloud interactions |
| Scheme specifics | Price and Rind (1992); Wong et al. (2013). Information from WRF-Chem. |
| External datasets | None, unless a dataset gives the IC to CG ratio |
| Information from MPAS-A | Land or ocean cover, level of neutral buoyancy (LNB), temperature, height, terrain height, grid cell area dA |
| Other inputs | Cloud-top height adjustment (km) |
| Information to MPAS-A | Intracloud flash rate IC_FR and cloud-to-ground flash rate CG_FR (flashes s-1) |
| Considerations | IC and CG can be split by a fixed ratio, by a dataset, or from storm properties. The coefficients A to E are prescribed. |
| Equations and references | total_FR = 3.44e-5 zkm^4.9 / 60 (land); 6.57e-6 zkm^4.9 / 60 (ocean); total_FR = total_FR dA / baseArea; IC_FR = total_FR (((A depth + B) depth C) depth + D) depth E; depth = (z_LNB - z_0C) 1e-3 + cldtop_adjustment, 5.5 < depth < 14 |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | discrete state update |
| Aerosol representation | none |
| Grid | 1-D (column) |
| Interdependence | needs convection scheme (LNB) |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `land_mask` | Land or ocean flag | 1 | to be defined | MPAS-A |
| `z_lnb` | Height of the level of neutral buoyancy | m | to be defined | MPAS-A |
| `temperature` | Air temperature profile | K | to be defined | MPAS-A |
| `height` | Height of each level | m | to be defined | MPAS-A |
| `terrain_height` | Terrain height | m | to be defined | MPAS-A |
| `cell_area` | Grid cell area | m2 | to be defined | MPAS-A |
| `cldtop_adjustment` | Cloud-top height adjustment | km | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `ic_flash_rate` | Intracloud flash rate | s-1 | to be defined | other parameterization (lightning NOx) |
| `cg_flash_rate` | Cloud-to-ground flash rate | s-1 | to be defined | other parameterization (lightning NOx) |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- Price and Rind (1992); Wong et al. (2013); information taken from WRF-Chem
