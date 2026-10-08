# biogenic/megan: MEGAN 3.0 biogenic VOC emissions

| | |
|---|---|
| **Status** | implemented |
| **Target class** | MEGAN Online Biogenic Emissions |
| **Contributors** | QUACS team (adapter); HuiW (upstream `MEGAN3_in_python`, MIT License, see [LICENSE](LICENSE)) |
| **Original source** | `MEGAN3_in_python`, a Python port of the site-scale MEGAN 3 model. The code in `src/` is a vendored and refactored copy. |

This scheme computes the emission flux of 19 biogenic VOC classes from vegetation for one grid cell and one time step.
A canopy model (`src/MEGCAN.py`) computes the sun and shade leaf temperature and light in each canopy layer.
The activity factors (`src/MEGVEA.py`) scale the emission factor of each class for light, temperature, leaf age, and optional stresses.
The adapter `biogenic_emission_megan_v3.py` gives one function call with scalar meteorology.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Biogenic VOC emissions, ozone chemistry, SOA |
| Scheme specifics | MEGAN 3.0 model |
| External datasets | Global maps of emission factors for some species. Climatological or past LAI, temperature, and solar radiation. In this example, `inputs/EF_LDF.csv` and `inputs/PFT_Fraction.csv` give the emission factors and the canopy fractions for one site. |
| Information from MPAS-A | LAI, surface solar radiation (as PPFD), 2-m temperature and humidity, wind speed, surface pressure, vegetation type (mapped to the 6 MEGAN canopy types), soil moisture, CO2 |
| Other inputs | VOC emission factors and light-dependent fractions (species and PFT dependent). Past-day statistics: 24-h mean temperature and PPFD, 10-day mean temperature, daily maximum and minimum temperature, daily maximum wind. |
| Information to MPAS-A | Emission flux of each VOC class, added to the concentration fields |
| Considerations | The model must map MEGAN classes to the CheMPAS-A species. The canopy model is part of this scheme; it can become a separate target parameterization. The host model must keep the past-day statistics. The upstream code averages the canopy state across PFTs before it computes the nonlinear activity responses. The adapter keeps this behavior. |
| Equations and references | Guenther et al. (2012, 2020); https://bai.ess.uci.edu/megan |

## Classification

| Dimension | Value |
|---|---|
| Complexity | complex |
| Solving strategy | supplies rates to the solver (emission flux) |
| Aerosol representation | none |
| Grid | 0-D (one surface cell, with an internal 1-D canopy) |
| Interdependence | standalone; needs past-day statistics from the host model |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
All inputs are keyword arguments of `compute_emissions` (`biogenic_emission_megan_v3`).
Scalar inputs are for one grid cell and one time step.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `day` | Day of year | day | scalar | MPAS-A |
| `hour` | Hour of day | h | scalar | MPAS-A |
| `latitude_deg` | Latitude | degrees | scalar | MPAS-A |
| `temperature_k` | Air temperature | K | scalar | MPAS-A |
| `pressure_pa` | Surface air pressure | Pa | scalar | MPAS-A |
| `wind_m_s` | Wind speed | m s-1 | scalar | MPAS-A |
| `relative_humidity_percent` | Relative humidity | % | scalar | MPAS-A |
| `ppfd` | Photosynthetic photon flux density above the canopy | umol m-2 s-1 | scalar | MPAS-A |
| `lai` | Leaf area index, this time step | m2 m-2 | scalar | MPAS-A |
| `previous_lai` | Leaf area index, previous period (leaf age response) | m2 m-2 | scalar | MPAS-A |
| `pft_fraction` | Fraction of each MEGAN canopy type (normalized inside) | 1 | `(6,)` | external dataset |
| `emission_factor` | Emission factor of each VOC class | nmol m-2 s-1 | `(19,)` | external dataset |
| `light_dependent_fraction` | Light-dependent fraction of each VOC class | 1 | `(19,)` | constant |
| `t24_k` | Mean temperature of the past 24 h | K | scalar | MPAS-A (history) |
| `p24` | Mean PPFD of the past 24 h | umol m-2 s-1 | scalar | MPAS-A (history) |
| `t10d_k` | Mean temperature of the past 10 days | K | scalar | MPAS-A (history) |
| `tmax_k` | Daily maximum temperature | K | scalar | MPAS-A (history) |
| `tmin_k` | Daily minimum temperature | K | scalar | MPAS-A (history) |
| `windmax_m_s` | Daily maximum wind speed | m s-1 | scalar | MPAS-A (history) |
| `kc_7d` | 7-day crop coefficient, for the soil moisture response | 1 | scalar | MPAS-A (history), optional |
| `swc30d` | Soil water content at 30 cm (checked, not used yet) | m3 m-3 | scalar | MPAS-A, optional |
| `settings` | `MeganSettings`: number of classes and layers, CO2, air quality index, solar constants | various | object | constant |
| `switches` | Dictionary that turns each activity response on or off | bool | dict | constant |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `flux` (return value) | Emission flux of each VOC class, in the order of `emission_factor` | nmol m-2 s-1 | `(19,)` | MPAS-A state (emission) or MICM rate parameter |

## Usage

```python
from quacs.biogenic.megan import MeganSettings, compute_emissions

flux = compute_emissions(day=180.0, hour=12.0, latitude_deg=45.4, temperature_k=298.15, ...)
```

| Example | What it does |
|---|---|
| `examples/driver_megan_v3.py` | Runs one grid cell for one time step with example meteorology and prints the flux of each VOC class |

## Tests

```bash
pytest quacs/biogenic/megan
```

## References

- Guenther, A. B., Jiang, X., Heald, C. L., Sakulyanontvittaya, T., Duhl, T., Emmons, L. K., and Wang, X. (2012). The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling biogenic emissions. Geosci. Model Dev., 5, 1471-1492. https://doi.org/10.5194/gmd-5-1471-2012
- Guenther, A., Jiang, X., Shah, T., Huang, L., Kemball-Cook, S., and Yarwood, G. (2020). Model of Emissions of Gases and Aerosol from Nature Version 3 (MEGAN3) for Estimating Biogenic Emissions. In: Air Pollution Modeling and its Application XXVI, Springer. (DOI to be added.)
- MEGAN home page: https://bai.ess.uci.edu/megan
