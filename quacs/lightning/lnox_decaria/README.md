# Lightning NOx for resolved convection

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Lightning-NOx Production for resolved convection |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Lightning-NOx effects on ozone, aerosols, aerosol-cloud interactions |
| Scheme specifics | DeCaria et al. (2000). Information from WRF-Chem. |
| External datasets | None |
| Information from MPAS-A | Land or ocean cover, temperature, air density, height, pressure, terrain height, 3-D radar reflectivity |
| Other inputs | IC_FR, CG_FR, N_IC, N_CG, temperatures of the peak IC and CG flashes, method to find cells where refl > refl_threshold |
| Information to MPAS-A | LNOx tendency (ppmv s-1) added to NO. Often also 1 or 2 LNOx tracers. |
| Considerations | NO must be in the MPAS scalar array. LNOx tracers, if used, must be in the scalar array and advected. |
| Equations and references | LNOx_IC_tend = IC_FR N_IC / cellcount fd_IC / B conv_factor where refl > refl_threshold (same form for CG). cellcount is the number of cells on each level with refl > refl_threshold; fd is a vertical distribution; B is a pressure factor. |

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
| `reflectivity` | 3-D radar reflectivity | dBZ | to be defined | MPAS-A |
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `air_density` | Air density | kg m-3 | to be defined | MPAS-A |
| `n_ic, n_cg` | Moles of NO for each IC and CG flash | mol | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `no_tendency` | NO tendency from lightning | ppmv s-1 | to be defined | MICM rate parameter |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- DeCaria et al. (2000), JGR; information taken from WRF-Chem
