# Reactive uptake of gases on mineral dust

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Heterogenous Chemistry |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Heterogeneous chemistry, tropospheric O3, dust uptake |
| Scheme specifics | Jacob (2000); Bauer et al. (2004); Dentener et al. (1996); Zhang et al. (2025); Tang et al. (2017) |
| External datasets | Uptake coefficients from laboratory measurements |
| Information from MPAS-A | Temperature, pressure, relative humidity |
| Other inputs | O3 mixing ratio, uptake coefficient gamma_O3 = 5e-5, dust effective radius and density, mean molecular speed of each gas |
| Information to MPAS-A | First-order loss rates for O3 (and other gases) as sink terms in the chemistry solver |
| Considerations | Also consider HNO3 (0.1), N2O5 (0.01-0.02), NO3 (1e-3), HO2 (0.2), and H2O2 (1e-3) uptake. Needs coupling to the gas-phase solver. |
| Equations and references | k = (1/4) gamma omega S_a for each gas |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | any |
| Grid | 0-D |
| Interdependence | needs dust surface area from the aerosol scheme |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `temperature` | Air temperature | K | to be defined | MPAS-A |
| `pressure` | Air pressure | Pa | to be defined | MPAS-A |
| `relative_humidity` | Relative humidity | 1 | to be defined | MPAS-A |
| `dust_surface_area` | Dust surface area density S_a | m2 m-3 | to be defined | other parameterization |
| `gamma` | Uptake coefficient for each gas | 1 | to be defined | constant |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `k_uptake` | First-order loss rate for each gas | s-1 | to be defined | MICM rate parameter |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- Jacob 2000, Bauer et al,, 2004, Dentener et al., 1996, Zhang et al., 2025, Tang et al., 2017. Reactive uptake of gases on mineral dust.
