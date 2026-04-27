#!/usr/bin/env python3
"""
Olson (2001) land-cover helper.

Loads resistance parameters from ``Olson_2001_Drydep_Inputs.nc`` and
exposes standalone functions for constructing :class:`LandCoverPatch` objects.

Import this module when you want to define land-cover patches from Olson
land-type indices rather than specifying resistance values manually::

    from quacs.drydep.lhs import olson_land_cover
    patches = [
        olson_land_cover.tropical_forest(fraction=0.6, lai=3.5),
        olson_land_cover.water(fraction=0.4),
    ]

Using :class:`LandCoverPatch` directly (without this module) is also valid
when resistance values are specified manually.
"""

from pathlib import Path

import numpy as np
import xarray as xr

from .box_model import DrydepCoefficients, LandCoverPatch

_data_dir = Path(__file__).parent.parent / "data"
_ds = xr.open_dataset(_data_dir / "Olson_2001_Drydep_Inputs.nc")

_IDEP = _ds.IDEP.values
_IRI = _ds.IRI.values.copy()
_IRI[2] = 200.0  # match current GEOS-Chem
_IRLU = _ds.IRLU.values
_IRAC = _ds.IRAC.values
_IRGSS = _ds.IRGSS.values
_IRGSO = _ds.IRGSO.values
_IRCLS = _ds.IRCLS.values
_IRCLO = _ds.IRCLO.values
_DRYCOEFF = _ds.DRYCOEFF.values
_IOLSON = _ds.IOLSON.values

# Pre-built coefficients instance
coefficients = DrydepCoefficients(drycoeff=_DRYCOEFF, iolson=_IOLSON)


def patch_from_olson(olson_id: int, fraction: float, lai: float) -> LandCoverPatch:
    """Build a :class:`LandCoverPatch` from an Olson 2001 land-type index.

    Parameters
    ----------
    olson_id : int
        Olson land-type index (0–72).
    fraction : float
        Fractional area within the grid cell (0–1).
    lai : float
        Leaf area index (m² m⁻²).
    """
    dep_id = int(_IDEP[olson_id])
    return LandCoverPatch(
        fraction=fraction,
        lai=lai,
        iri=float(_IRI[dep_id - 1]),
        irlu=float(_IRLU[dep_id - 1]),
        irac=float(_IRAC[dep_id - 1]),
        irgss=float(_IRGSS[dep_id - 1]),
        irgso=float(_IRGSO[dep_id - 1]),
        ircls=float(_IRCLS[dep_id - 1]),
        irclo=float(_IRCLO[dep_id - 1]),
    )


def tropical_forest(fraction: float = 1.0, lai: float = 0.71) -> LandCoverPatch:
    """Tropical/subtropical broadleaf forest (Olson type 29)."""
    return patch_from_olson(29, fraction, lai)


def boreal_forest(fraction: float = 1.0, lai: float = 0.61) -> LandCoverPatch:
    """Boreal/temperate needleleaf forest (Olson type 33)."""
    return patch_from_olson(33, fraction, lai)


def savanna(fraction: float = 1.0, lai: float = 0.74) -> LandCoverPatch:
    """Savanna / shrubland (Olson type 43)."""
    return patch_from_olson(43, fraction, lai)


def water(fraction: float = 1.0, lai: float = 0.1) -> LandCoverPatch:
    """Water / ocean surface (Olson type 0)."""
    return patch_from_olson(0, fraction, lai)


def cropland(fraction: float = 1.0, lai: float = 0.5) -> LandCoverPatch:
    """Cropland / agriculture (Olson type 9)."""
    return patch_from_olson(9, fraction, lai)
