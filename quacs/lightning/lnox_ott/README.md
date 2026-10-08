# Lightning NOx for parameterized convection

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Lightning-NOx Production for parameterized convection |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Lightning-NOx effects on ozone, aerosols, aerosol-cloud interactions |
| Scheme specifics | Ott et al. (2010). Information from WRF-Chem. |
| External datasets | None |
| Information from MPAS-A | Land or ocean cover, latitude, air density, height, terrain height |
| Other inputs | IC_FR, CG_FR, moles of NO for each IC flash (N_IC) and each CG flash (N_CG) |
| Information to MPAS-A | LNOx tendency (ppmv s-1) added to NO. Often also 1 or 2 LNOx tracers. |
| Considerations | NO must be in the MPAS scalar array. LNOx tracers, if used, must be in the scalar array and advected. |
| Equations and references | LNOx_tend = (IC_FR N_IC + CG_FR N_CG) vertdist(lat, land) mass_air. There are 4 prescribed vertical distributions. |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | none |
| Grid | 1-D (column) |
| Interdependence | needs a flash rate parameterization |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `ic_flash_rate` | Intracloud flash rate | s-1 | to be defined | other parameterization |
| `cg_flash_rate` | Cloud-to-ground flash rate | s-1 | to be defined | other parameterization |
| `land_mask` | Land or ocean flag | 1 | to be defined | MPAS-A |
| `latitude` | Latitude | degrees | to be defined | MPAS-A |
| `air_density` | Air density | kg m-3 | to be defined | MPAS-A |
| `height` | Height of each level | m | to be defined | MPAS-A |
| `n_ic, n_cg` | Moles of NO for each IC and CG flash | mol | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `no_tendency` | NO tendency from lightning | ppmv s-1 | to be defined | MICM rate parameter |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- Ott et al. (2010), JGR; information taken from WRF-Chem
