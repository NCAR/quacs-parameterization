# LS99-FENGSHA dust emission

| | |
|---|---|
| **Status** | placeholder |
| **Target class** | Dust Emission Scheme |
| **Contributors** | Open. Add your name when you start this scheme. |
| **Original source** | To be added |

This folder is a placeholder. The content below comes from the QUACS target parameterization list.
The inputs and outputs are a first draft. Change them when you implement the scheme.

## Target parameterization template

| Field | Value |
|---|---|
| Research area | Emissions, aerosols, dust storms, radiative impacts |
| Scheme specifics | LS99-FENGSHA dust emission scheme (Ma et al., 2019; Foroutan et al., 2017; CMAQ 5.5) |
| External datasets | Soil texture, MODIS FPAR, land-use classification. A lookup table (Kang et al., 2011) gives the fine-particle fraction f, plastic pressure p, bulk soil density, particle density, and the sandblasting parameters C_alpha and C_beta. |
| Information from MPAS-A | u10m, volumetric soil moisture of the top layer, snow cover, FPAR |
| Other inputs | A_n = 0.0123, Gamma = 1.65e-4 kg s-2. Soil moisture correction from Fécan et al. (1999). Roughness length z0 from FPAR (Foroutan et al., 2017). |
| Information to MPAS-A | Vertical dust flux F [x, y], split into accumulation mode (mass fraction 0.07, Dp,g = 1.39 um) and coarse mode (mass fraction 0.93, Dp,g = 5.26 um), sigma_g = 2 |
| Considerations | The scheme computes u* and z0 itself. They do not come from MPAS-A. Dust has dry deposition. |
| Equations and references | u_t0(D) = sqrt(A_n rho_p g D / rho_air + Gamma / (rho_air D)); u*_t = u_t0 sqrt(1 + 1.21 (w - w')^0.68) for w > w'; z0 = 0.96 lambda^1.07 h_veg; u* = kappa u10 / ln(z_ref / z0); Q = (rho_air / g) u*^3 (1 + u*_t/u*) (1 - u*_t^2/u*^2); F = alpha Q |

## Classification

| Dimension | Value |
|---|---|
| Complexity | simple |
| Solving strategy | discrete state update |
| Aerosol representation | modal |
| Grid | 0-D |
| Interdependence | standalone |

## Inputs

| Name | Description | Units | Shape | Source |
|---|---|---|---|---|
| `u10m` | 10-m wind speed | m s-1 | to be defined | MPAS-A |
| `soil_moisture` | Volumetric soil moisture, top layer | m3 m-3 | to be defined | MPAS-A |
| `snow_cover` | Snow cover fraction | 1 | to be defined | MPAS-A |
| `fpar` | Fraction of absorbed photosynthetically active radiation | 1 | to be defined | MPAS-A |
| `soil_texture` | Soil texture class | 1 | to be defined | external dataset |
| `land_use` | Land-use class | 1 | to be defined | external dataset |

## Outputs

| Name | Description | Units | Shape | Destination |
|---|---|---|---|---|
| `dust_flux_accumulation` | Vertical dust flux, accumulation mode | kg m-2 s-1 | to be defined | MPAS-A state |
| `dust_flux_coarse` | Vertical dust flux, coarse mode | kg m-2 s-1 | to be defined | MPAS-A state |

## References

This section shows each reference exactly as its source gives it. Nobody added titles, DOIs, or page numbers. Check them before you cite this scheme.

From the target parameterization list:

- Ma et al., 2019
- Foroutan et al., 2017 (also the roughness length z0 from FPAR)
- CMAQ 5.5
- Kang et al., 2011 (lookup table of soil parameters)
- Fecan et al., 1999 (soil moisture correction)
