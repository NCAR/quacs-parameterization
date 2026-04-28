"""
LHS parameter space definition for dry deposition sensitivity analysis.

Defines the meteorological and land-cover dimensions sampled by the Latin
Hypercube driver, and provides ``build_patches_from_sample`` to convert raw
LHS samples into ``LandCoverPatch`` lists.

Keeping this separate from the driver script lets tests import these
definitions without pulling in matplotlib, argparse, or musica.
"""

import numpy as np

from quacs.drydep.lhs import olson_land_cover

# Each entry: (met_key, lower_bound, upper_bound, description)
MET_DIMS = [
    ("TC0", 270.0, 310.0, "Temperature (K)"),
    ("PRESSU", 5.0e4, 1.05e5, "Pressure (Pa)"),
    ("USTAR", 0.05, 1.0, "Friction velocity (m s⁻¹)"),
    ("RADIAT", 0.0, 900.0, "Shortwave radiation (W m⁻²)"),
    ("CFRAC", 0.0, 1.0, "Cloud fraction"),
    ("SUNCOS_MID", 0.01, 1.0, "cos(solar zenith angle)"),
    ("AZO", 0.01, 5.0, "Roughness height (m)"),
    ("AIRDEN", 0.5, 1.3, "Air density (kg m⁻³)"),
    ("HFLUX", -50.0, 100.0, "Sensible heat flux (W m⁻²)"),
    ("ALBD", 0.05, 0.9, "Surface albedo"),
    ("BXHEIGHT", 50.0, 2000.0, "Boundary layer height (m)"),
    ("U10M", -10.0, 10.0, "Zonal 10-m wind (m s⁻¹)"),
    ("V10M", -10.0, 10.0, "Meridional 10-m wind (m s⁻¹)"),
]

# Each entry: (olson_land_cover function name, lai_min, lai_max)
LC_ARCHETYPES = [
    ("tropical_forest", 0.1, 7.0),
    ("boreal_forest", 0.1, 5.0),
    ("savanna", 0.1, 3.0),
    ("cropland", 0.1, 4.0),
    ("water", 0.0, 0.0),  # LAI fixed at 0 for water
]
N_LC_ARCHETYPES = len(LC_ARCHETYPES)

_WATER_IDX = next(i for i, lc in enumerate(LC_ARCHETYPES) if lc[0] == "water")


def build_patches_from_sample(frac_raw: np.ndarray, lai_vals: np.ndarray):
    """Convert raw LHS fractions and LAI values into a list of LandCoverPatch.

    Parameters
    ----------
    frac_raw : ndarray, shape (N_LC_ARCHETYPES,)
        Un-normalised fractions; normalised internally to sum to 1.
    lai_vals : ndarray, shape (N_LC_ARCHETYPES,)
        Leaf area index for each archetype (water LAI is forced to 0).
    """
    fracs = frac_raw / frac_raw.sum()
    patches = []
    for i, (name, _lai_min, _lai_max) in enumerate(LC_ARCHETYPES):
        lai = 0.0 if name == "water" else float(lai_vals[i])
        constructor = getattr(olson_land_cover, name)
        patches.append(constructor(fraction=float(fracs[i]), lai=lai))
    return patches
