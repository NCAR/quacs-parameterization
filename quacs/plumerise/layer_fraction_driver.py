"""Step 3 driver: build `layer_fraction[k]`."""

from __future__ import annotations

import numpy as np

from quacs.plumerise.tables import lookup_vegetation_split

PYROCB_FRACTION = 0.5


def resolve_mass_split(
    vegetation_class,
    pyrocb_flag,
):
    split, _source_key = lookup_vegetation_split(vegetation_class)
    remaining_fraction = 1.0
    pyrocb_fraction = 0.0
    if pyrocb_flag:
        pyrocb_fraction = PYROCB_FRACTION
        remaining_fraction = 1.0 - pyrocb_fraction
    smoldering_weight = remaining_fraction * split.smoldering_fraction
    flaming_weight = remaining_fraction * split.flaming_fraction
    pyrocb_weight = pyrocb_fraction
    return smoldering_weight, flaming_weight, pyrocb_weight


def zero_layer_fraction(layer_bounds):
    return np.zeros(len(layer_bounds) - 1, dtype=float)


def top_hat_overlap_fraction(
    layer_bounds,
    height_base_m,
    height_top_m,
):
    if len(layer_bounds) < 2:
        raise ValueError("layer_bounds must contain at least two bounds")
    if np.any(np.diff(layer_bounds) <= 0.0):
        raise ValueError("layer_bounds must be strictly increasing")
    if height_top_m <= height_base_m:
        raise ValueError("height_top_m must exceed height_base_m")

    layer_low = layer_bounds[:-1]
    layer_high = layer_bounds[1:]
    overlaps = np.maximum(
        0.0,
        np.minimum(layer_high, height_top_m) - np.maximum(layer_low, height_base_m),
    )
    total = float(np.sum(overlaps))
    if total <= 0.0:
        return zero_layer_fraction(layer_bounds)
    return overlaps / total


def surface_layer_fraction(layer_bounds):
    fraction = zero_layer_fraction(layer_bounds)
    if fraction.size:
        fraction[0] = 1.0
    return fraction


def combine_layer_fractions(
    smoldering_weight,
    flaming_weight,
    pyrocb_weight,
    smoldering_layer_fraction,
    flaming_layer_fraction,
    pyrocb_layer_fraction,
):
    if not (
        len(smoldering_layer_fraction)
        == len(flaming_layer_fraction)
        == len(pyrocb_layer_fraction)
    ):
        raise ValueError("layer fractions must share length")
    return (
        smoldering_weight * smoldering_layer_fraction
        + flaming_weight * flaming_layer_fraction
        + pyrocb_weight * pyrocb_layer_fraction
    )


def layer_fraction_driver(
    z,
    vegetation_class,
    injectH_base_m,
    injectH_top_m,
    injectH_pyroCb_base_m,
    injectH_pyroCb_top_m,
    pyrocb_flag,
):
    """Return direct `layer_fraction[k]` from height and split variables."""

    if len(z) < 2:
        raise ValueError("z must contain at least two levels")
    dz_top = z[-1] - z[-2]
    layer_bounds = np.concatenate([z, np.array([z[-1] + dz_top])])
    smoldering_weight, flaming_weight, pyrocb_weight = resolve_mass_split(
        vegetation_class,
        pyrocb_flag,
    )
    smoldering_layer_fraction = surface_layer_fraction(layer_bounds)
    flaming_layer_fraction = top_hat_overlap_fraction(
        layer_bounds,
        injectH_base_m,
        injectH_top_m,
    )
    pyrocb_layer_fraction = zero_layer_fraction(layer_bounds)
    if (
        pyrocb_flag
        and injectH_pyroCb_base_m is not None
        and injectH_pyroCb_top_m is not None
    ):
        pyrocb_layer_fraction = top_hat_overlap_fraction(
            layer_bounds,
            injectH_pyroCb_base_m,
            injectH_pyroCb_top_m,
        )
    layer_fraction = combine_layer_fractions(
        smoldering_weight,
        flaming_weight,
        pyrocb_weight,
        smoldering_layer_fraction,
        flaming_layer_fraction,
        pyrocb_layer_fraction,
    )
    return layer_fraction
