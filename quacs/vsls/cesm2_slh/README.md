# vsls/cesm2_slh: sea-salt bromine recycling (CESM2-SLH)

| | |
|---|---|
| **Status** | in progress |
| **Target class** | Unique Reaction Rate - Halogenated very short-lived substances (VSLS) |
| **Contributors** | Saeideh Mohammadi (University of Iowa) |
| **Original source** | CAM `src/chemistry/mozart/mo_usrrxt.F90` ([ESCOMP/CAM](https://github.com/ESCOMP/CAM), `cam_development` branch) |

This scheme is the first component of the VSLS target parameterization.
It computes the pseudo-first-order rate constants for the uptake of reactive bromine on sea-salt aerosol.
The uptake releases Br2 and BrCl to the gas phase (sea-salt dehalogenation), as in CESM2-SLH.

The Python code supplies only the non-standard rate constants. In CAM, these rate constants are in `mo_usrrxt.F90`.
The MICM and TUV-x configuration files give the species, the standard reactions, and the photolysis.
MICM does the time integration.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Tropospheric and lower-stratospheric ozone, halogens, methane lifetime, oxidation capacity |
| Scheme specifics | CESM2-SLH sea-salt dehalogenation, CAM reactions `het_ss_0`, `het_ss_1`, `het_ss_2` (`mo_usrrxt.F90`). Part of the VSLS scheme (Ordóñez et al., 2012; Fernandez et al., 2014). |
| External datasets | None for this component |
| Information from MPAS-A | Air temperature, mid-layer pressure, latitude, ocean plus sea-ice fraction, calendar day, sea-salt surface area density |
| Other inputs | Uptake coefficients (gamma) and product yields, depletion factor bounds, `SSAdehal_ScalingFactor` (CAM namelist, default 1.0) |
| Information to MPAS-A | First-order rate constants (s-1) that go to MICM as `USER.het_ss_0`, `USER.het_ss_1`, `USER.het_ss_2`. MICM returns updated BrONO2, BrNO2, HOBr, Br2, and BrCl. |
| Considerations | See [Considerations](#considerations) |
| Equations and references | See [Equations](#equations) and [References](#references) |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | supplies rates to the solver |
| Aerosol representation | any (needs only the sea-salt surface area density) |
| Grid | 0-D |
| Interdependence | needs sea-salt surface area density from an aerosol scheme |

## Inputs

Source is one of: MPAS-A, external dataset, constant, other parameterization.
All arrays can be scalars or arrays that broadcast together.

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `temperature_k` | Air temperature | K | scalar or `(ncell,)` | MPAS-A |
| `pressure_pa` | Mid-layer pressure | Pa | scalar or `(ncell,)` | MPAS-A |
| `latitude_deg` | Latitude | degrees north | scalar or `(ncell,)` | MPAS-A |
| `calday` | Calendar day of the year (CAM convention) | day | scalar or `(ncell,)` | MPAS-A |
| `sad_seasalt_m2_m3` | Sea-salt surface area density | m2 m-3 | scalar or `(ncell,)` | other parameterization |
| `ocean_ice_fraction` | Ocean plus sea-ice fraction of the grid cell | 1 | scalar or `(ncell,)` | MPAS-A |
| `scaling_factor` | CAM `SSAdehal_ScalingFactor` | 1 | scalar | constant (default 1.0) |

## Outputs

Destination is one of: MPAS-A state, MICM rate parameter, diagnostic.

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `het_ss_0` | Rate constant for BrONO2 -> 0.65 Br2 + 0.35 BrCl | s-1 | same as the inputs | MICM rate parameter `USER.het_ss_0` |
| `het_ss_1` | Rate constant for BrNO2 -> 0.65 Br2 + 0.35 BrCl | s-1 | same as the inputs | MICM rate parameter `USER.het_ss_1` |
| `het_ss_2` | Rate constant for HOBr -> 0.65 Br2 + 0.35 BrCl | s-1 | same as the inputs | MICM rate parameter `USER.het_ss_2` |

`seasalt_rate_constants()` returns these outputs as a dictionary.
`build_mechanism()` returns the MICM mechanism that uses them.

## Reactions

| CAM tag | Reaction | gamma |
|---|---|---|
| `het_ss_0` | BrONO2 -> 0.65 Br2 + 0.35 BrCl | 0.010 |
| `het_ss_1` | BrNO2 -> 0.65 Br2 + 0.35 BrCl | 0.005 |
| `het_ss_2` | HOBr -> 0.65 Br2 + 0.35 BrCl | 0.0125 |

The comments in the CAM code tell that these gamma values are half of the original CESM2 values.

## Equations

The rate constant uses free-molecular uptake, as in CAM (units s-1):

    k = S * 0.25 * gamma * v_mean * A_ss * DF * mask

    v_mean = sqrt(8 * R * T / (pi * M))      mean molecular speed (m/s)

A_ss is the sea-salt surface area density (m2/m3). DF is the bromine depletion factor.
S is `SSAdehal_ScalingFactor`. The mask is 0 or 1.

The depletion factor is:

    DF = 0.5                                                                north of 30S
    DF = 0.9 + 0.5 * (0.3 - 0.9) * (sin((calday/182.5 - 0.5) * pi) + 1)    south of 30S

The mask is 0 at pure land points north of 60S, because CAM sets the sea-salt SAD to zero there.
The mask is also 0 at pressures below 300 hPa. At all other points, the mask is 1.

## Considerations

- The VSLS template on the target list does not give all the inputs that this component needs:
  - Sea-salt SAD. CheMPAS-A v1 has no aerosols. GOCART-2G carries sea salt, but it links only in v2. The source of sea-salt SAD in CheMPAS-A is an open question.
  - Pressure, for the 300 hPa cutoff.
  - Calendar day, for the depletion factor south of 30S.
- These 3 reactions do not need chemical concentrations. The ice and sulfate reactions will need them.
- Gas-phase bromine is not conserved, on purpose. Each Br atom that the sea salt takes up releases 1.65 Br atoms to the gas phase. The extra 0.65 Br comes from bromide in the sea salt, which is not a tracked species. A test checks this.
- CAM uses SAD in cm2/cm3. 1 cm2/cm3 = 100 m2/m3.
- The CAM free-molecular form does not include the gas-phase diffusion limit. The MICM `Surface` reaction includes this limit, but it has no depletion factor, mask, or scaling factor. Thus this scheme uses `UserDefined` reactions, to stay consistent with CAM.
- CAM hard-codes the speed prefactor 100*sqrt(8R/(pi*M)), for example 1.47e3 for HOBr. This code calculates it from the molar mass. The two values agree within 1%, and a test checks this.
- The heterogeneous chemistry target ([`heterogeneous_chem/dust_uptake`](../../heterogeneous_chem/dust_uptake/)) uses the same form k = 0.25 * gamma * v * A. The two schemes can share one uptake function.

### Not implemented yet (open questions)

- HOBr + HCl and HOBr + HBr on ice and sulfate. These reactions consume a second gas. CAM uses a limiting-reactant rate that divides by the HOBr or HCl number density, so the rate constant depends on the concentrations. This needs a scripted bimolecular rate. The MICM `Surface` reaction cannot do it.
- pH-dependent and halide-dependent gamma with variable Br2/BrCl branching, as in GEOS-Chem.
- The rest of the VSLS scheme: the mechanism configuration, the emissions, the deposition of the new species, and the new TUV-x photolysis rates.

## Usage

```python
from quacs.vsls.cesm2_slh import build_mechanism, rate_parameter_key, seasalt_rate_constants

rates = seasalt_rate_constants(
    temperature_k=288.0,
    pressure_pa=1.0e5,
    latitude_deg=40.0,
    calday=180.0,
    sad_seasalt_m2_m3=5.0e-5,
    ocean_ice_fraction=1.0,
)
mechanism, species = build_mechanism()
```

| Example | What it does |
|---|---|
| `examples/musica_seasalt_box_model.py` | Runs the 3 reactions in one marine cell for 24 h with MICM and compares the result with the analytic solution |
| `examples/run_vsls_seasalt.ipynb` | Plots the rate constants and the MICM box model result |

## Tests

```bash
pytest quacs/vsls/cesm2_slh
```

## References

CAM source files (ESCOMP/CAM, `cam_development` branch):

- `src/chemistry/mozart/mo_usrrxt.F90`: gamma values, DF, masks, rate form
- `src/chemistry/mozart/mo_slh_routines.F90`: default `SSAdehal_ScalingFactor`
- `src/chemistry/pp_trop_strat_mam4_slh/chem_mech.in`: products and yields

Literature that the QUACS VSLS template cites:

- Ordóñez, C., et al. (2012). Bromine and iodine chemistry in a global chemistry-climate model: description and evaluation of very short-lived oceanic sources. Atmos. Chem. Phys., 12, 1423-1447. https://doi.org/10.5194/acp-12-1423-2012
- Fernandez, R. P., et al. (2014). Bromine partitioning in the tropical tropopause layer: implications for stratospheric injection. Atmos. Chem. Phys., 14, 13391-13410. https://doi.org/10.5194/acp-14-13391-2014
- Yang, X., et al. (2005). Tropospheric bromine chemistry and its impacts on ozone: A model study. J. Geophys. Res., 110, D23311. https://doi.org/10.1029/2005JD006244
