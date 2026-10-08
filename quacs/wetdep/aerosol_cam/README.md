# CAM aerosol wet deposition (wetdepa_v2)

| | |
|---|---|
| **Status** | implemented |
| **Target class** | Resolved scale Aerosol Wet Deposition |
| **Contributors** | Sue Park |
| **Original source** | CAM `src/chemistry/aerosol/wetdep.F90`, subroutines `clddiag` and `wetdepa_v2` (https://github.com/ESCOMP/CAM) |

This scheme is a Python port of the CAM (CESM) bulk aerosol wet scavenging routines.
It computes the removal of an aerosol tracer by in-cloud (nucleation) scavenging and by below-cloud (impaction) scavenging.
It releases the scavenged mass again in the layers where the precipitation evaporates.
The port covers the non-cloudborne branch of `wetdepa_v2` only.

The target list names the PNNL MOSAIC sectional scheme for this class (see [`aerosol_mosaic`](../aerosol_mosaic/)).
This CAM bulk scheme is a second scheme in the same class.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Aerosol life cycle, air quality, climate |
| Scheme specifics | CAM `wetdepa_v2`, bulk (non-modal) aerosol, non-cloudborne branch. In-cloud scavenging for stratiform and convective clouds. Below-cloud scavenging with the Dana and Hales (1976) coefficient. Resuspension where precipitation evaporates. |
| External datasets | None |
| Information from MPAS-A | Layer pressure thickness, total and convective cloud fraction, stratiform cloud water, in-cloud convective water, stratiform and convective precipitation production and evaporation, detrained convective water, tracer mass mixing ratio, time step |
| Other inputs | Solubility factors (`sol_fact`, `sol_facti`, `sol_factic`), below-cloud scavenging coefficient `scavcoef` (default 0.1 mm-1) |
| Information to MPAS-A | Tracer tendency, surface wet deposition flux. For MICM: first-order loss rates `k_strat + k_conv` and a resuspension source for each level. |
| Considerations | The scheme needs the vertical structure of the column. It sweeps from the model top (k = 0) to the surface, because the precipitation flux into each layer comes from the layers above. The cloud-borne branch, the MAM size-dependent below-cloud coefficients, and the ice and snow factors are not ported. |
| Equations and references | Rate equations are in the docstring of `wetdep_cam.py`. Rasch et al. (2000); Barth et al. (2000); Dana and Hales (1976). |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver (first-order loss and resuspension source); `wetdepa_v2` also gives the CAM explicit tendency |
| Aerosol representation | bulk |
| Grid | 1-D (column) |
| Interdependence | needs cloud and precipitation fields from the host microphysics and convection schemes |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
The array layout follows CAM: `[ncol, nlev]`. Level `k = 0` is the model top and `k = nlev - 1` is next to the surface.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `q` | Tracer mass mixing ratio | kg kg-1 | `(ncol, nlev)` | MPAS-A |
| `pdel` | Layer pressure thickness | Pa | `(ncol, nlev)` | MPAS-A |
| `cldt` | Total cloud fraction | 1 | `(ncol, nlev)` | MPAS-A |
| `cldc` | Convective cloud fraction | 1 | `(ncol, nlev)` | MPAS-A |
| `cwat` | Stratiform cloud water | kg kg-1 | `(ncol, nlev)` | MPAS-A |
| `precs` | Stratiform precipitation production rate | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A |
| `evaps` | Stratiform precipitation evaporation rate | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A |
| `conicw` | In-cloud convective water | kg kg-1 | `(ncol, nlev)` | MPAS-A |
| `cmfdqr` | Convective precipitation production rate | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A |
| `evapc` | Convective precipitation evaporation rate | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A |
| `dlf` | Detrained convective water (optional, default 0) | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A |
| `deltat` | Time step | s | scalar | MPAS-A |
| `sol_fact` | Solubility factor: fraction of the tracer that scavenging can remove | 1 | scalar | constant |
| `sol_facti`, `sol_factic` | Solubility factors for stratiform and convective in-cloud scavenging (optional, default `sol_fact`) | 1 | scalar | constant |
| `scavcoef` | Below-cloud scavenging coefficient (default `SCAVCOEF_DEFAULT` = 0.1) | mm-1 | scalar | constant |
| `cldvst`, `cldvcu` | Precipitating cloud area, stratiform and convective (optional; `cloud_volume_diag` computes them if not given) | 1 | `(ncol, nlev)` | other parameterization |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

`wetdepa_v2` returns a dictionary with these keys:

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `scavt` | Total tracer tendency, with resuspension | kg kg-1 s-1 | `(ncol, nlev)` | MPAS-A state |
| `ic_scavt` | In-cloud removal part of the tendency (not positive) | kg kg-1 s-1 | `(ncol, nlev)` | diagnostic |
| `bc_scavt` | Below-cloud removal part of the tendency (not positive) | kg kg-1 s-1 | `(ncol, nlev)` | diagnostic |
| `resusp` | Resuspension source (not negative) | kg kg-1 s-1 | `(ncol, nlev)` | diagnostic |
| `sfc_flux` | Surface wet deposition flux | kg m-2 s-1 | `(ncol,)` | diagnostic |
| `coef` | Output of `column_scavenging_coefficients` (below) | various | dict | diagnostic |

`column_scavenging_coefficients` returns the rates that do not depend on `q`. The MUSICA examples give these to MICM:

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `k_strat`, `k_conv` | First-order removal rates after the limiter, stratiform and convective | s-1 | `(ncol, nlev)` | MICM rate parameter (`LOSS.*`) |
| `fracev`, `fracev_cu` | Fraction of the precipitation that evaporates in the layer | 1 | `(ncol, nlev)` | diagnostic |
| `precabs`, `precabc` | Precipitation flux into the layer from above | kg m-2 s-1 | `(ncol, nlev)` | diagnostic |
| `cldvst`, `cldvcu` | Precipitating cloud area | 1 | `(ncol, nlev)` | diagnostic |
| `pdog` | Layer mass, `pdel / g` | kg m-2 | `(ncol, nlev)` | diagnostic |

`column_step_exponential(q, coef, deltat)` returns `q_new` (kg kg-1), the resuspension source `s` (kg kg-1 s-1, given to MICM as `EMIS.*`), and the step-mean `sfc_flux` (kg m-2 s-1).

## Usage

```python
from quacs.wetdep.aerosol_cam import wetdepa_v2

result = wetdepa_v2(
    q,
    pdel,
    cldt,
    cldc,
    cwat,
    precs,
    evaps,
    conicw,
    cmfdqr,
    evapc,
    deltat=1800.0,
    sol_fact=0.3,
)
q_new = q + result["scavt"] * 1800.0  # CAM explicit update
```

| Example | What it does |
|---|---|
| `examples/single_column_demo.py` | Prints the rates and the tendency for one column with one stratiform cloud layer |
| `examples/musica_box_model.py` | Gives the first-order loss rate of one layer to a MICM box model and compares MICM with the exact and the CAM explicit solutions |
| `examples/musica_column_model.py` | Runs each level as one MICM grid cell, with loss and resuspension, and checks the column mass budget |

## Tests

```bash
pytest quacs/wetdep/aerosol_cam
```

The tests check zero removal with no precipitation, the signs of the terms, linearity in `q`, the column mass budget, full resuspension, the precipitating cloud area, hand calculations of the in-cloud and below-cloud rates, the limiter, and the agreement of the exact step with the CAM explicit update for a small time step.

## References

- Rasch, P. J., Feichter, J., Law, K., et al. (2000). A comparison of scavenging and deposition processes in global models: results from the WCRP Cambridge Workshop of 1995. Tellus B, 52, 1025-1056. https://doi.org/10.3402/tellusb.v52i4.17085
- Barth, M. C., Rasch, P. J., Kiehl, J. T., Benkovitz, C. M., and Schwartz, S. E. (2000). Sulfur chemistry in the National Center for Atmospheric Research Community Climate Model: Description, evaluation, features, and sensitivity to aqueous chemistry. J. Geophys. Res., 105, 1387-1415. https://doi.org/10.1029/1999JD900773
- Dana, M. T. and Hales, J. M. (1976). Statistical aspects of the washout of polydisperse aerosols. Atmos. Environ., 10, 45-50. https://doi.org/10.1016/0004-6981(76)90258-1
