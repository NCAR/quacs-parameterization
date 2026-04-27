#!/usr/bin/env python3
"""
Dry deposition box model helper.

Provides:
  - ``LandCoverPatch``   – dataclass carrying per-patch land cover parameters
  - ``DrydepCoefficients`` – dataclass for Wesely scheme algorithm coefficients
  - ``compute_drydep_rate`` – vectorized; accepts a list of
                              ``musica.mechanism_configuration.Species``, optional
                              scalar or array met conditions, and a list of
                              ``LandCoverPatch``.  Returns a dict mapping species
                              name → ndarray of first-order loss rates, or a tuple
                              in single-species scalar mode for backward compat.

Each ``mc.Species`` must carry the following keys in ``other_properties``:
  ``"henrys_law_constant"`` – effective Henry's law constant at pH 7 (M atm⁻¹)
  ``"reactivity"``          – biological reactivity factor F0 (dimensionless, 0–1)

These values are stored as strings because the MUSICA C++ binding requires
``Mapping[str, str]``; use ``float(sp.other_properties[key])`` to recover floats.

Standard species presets (``hg0``, ``so2``, ``o3``) are defined in
``musica_box_model``.  Olson 2001 land-cover helpers are in
``olson_land_cover``; this module has no knowledge of specific species or
land-cover data sources.
"""

from dataclasses import dataclass
from typing import List, Union

import musica.mechanism_configuration as mc
import numpy as np

from quacs.drydep.simple.drydep_functions import DEPVEL, METERO


# ── DrydepCoefficients ───────────────────────────────────────────────────────

@dataclass
class DrydepCoefficients:
    """Algorithm coefficients required by the DEPVEL dry deposition routine.

    These are properties of the resistance parameterisation scheme, not of
    any particular land cover.  Use ``olson_land_cover.coefficients`` to
    obtain the standard Olson 2001 values.

    Parameters
    ----------
    drycoeff : np.ndarray
        Wesely (1989) dimensionless reactivity/solubility coefficients.
    iolson : np.ndarray
        Olson land-type to dep-type mapping array used by DEPVEL.
    """

    drycoeff: np.ndarray
    iolson: np.ndarray


# ── LandCoverPatch ────────────────────────────────────────────────────────────

@dataclass
class LandCoverPatch:
    """Parameters describing one land cover type within a grid cell.

    A grid cell is represented as a list of ``LandCoverPatch`` objects whose
    ``fraction`` values sum to 1.

    Parameters
    ----------
    fraction : float
        Fractional area of this patch in the grid cell (0–1).
    lai : float
        Leaf area index (m² m⁻²).
    iri : float
        Minimum stomatal resistance (s m⁻¹).
    irlu : float
        Cuticular resistance (s m⁻¹, per unit LAI).
    irac : float
        Canopy aerodynamic resistance (s m⁻¹).
    irgss : float
        Soil/ground resistance for SO₂ analogue (s m⁻¹).
    irgso : float
        Soil/ground resistance for O₃ analogue (s m⁻¹).
    ircls : float
        Lower-canopy resistance for SO₂ analogue (s m⁻¹).
    irclo : float
        Lower-canopy resistance for O₃ analogue (s m⁻¹).
    """

    fraction: float
    lai: float
    iri: float
    irlu: float
    irac: float
    irgss: float
    irgso: float
    ircls: float
    irclo: float


# ── Internal helper ──────────────────────────────────────────────────────────

def _patches_to_depvel_arrays(patches: List[LandCoverPatch]):
    """Convert a list of LandCoverPatch objects to DEPVEL-compatible arrays."""
    n = len(patches)
    IREG = n
    ILAND = np.zeros(73, dtype=int)
    IUSE = np.zeros(73)
    XLAI = np.zeros(73)
    idep_arr = np.zeros(73, dtype=int)
    iri_arr = np.zeros(11)
    irlu_arr = np.zeros(11)
    irac_arr = np.zeros(11)
    irgss_arr = np.zeros(11)
    irgso_arr = np.zeros(11)
    ircls_arr = np.zeros(11)
    irclo_arr = np.zeros(11)

    for i, patch in enumerate(patches):
        ILAND[i] = i
        IUSE[i] = patch.fraction * 1000.0
        XLAI[i] = patch.lai
        idep_arr[i] = i + 1
        iri_arr[i] = patch.iri
        irlu_arr[i] = patch.irlu
        irac_arr[i] = patch.irac
        irgss_arr[i] = patch.irgss
        irgso_arr[i] = patch.irgso
        ircls_arr[i] = patch.ircls
        irclo_arr[i] = patch.irclo

    return (
        IREG, ILAND, IUSE, XLAI,
        idep_arr, iri_arr, irlu_arr, irac_arr,
        irgss_arr, irgso_arr, ircls_arr, irclo_arr,
    )


# ── compute_drydep_rate ───────────────────────────────────────────────────────

def compute_drydep_rate(
    species: List[mc.Species],
    *,
    box_height_m: float = 1.0,
    met: dict,
    land_cover: Union[List[LandCoverPatch], List[List[LandCoverPatch]]],
    coefficients: DrydepCoefficients,
):
    """Compute dry deposition velocities and first-order loss rates.

    **Single-species scalar mode** (one species, scalar met)::

        dvel_cms, k_s = compute_drydep_rate([hg0], met={...}, land_cover=[...],
                                             coefficients=olson_land_cover.coefficients)

    **Vectorized mode** (N cells or multiple species)::

        result = compute_drydep_rate(
            [hg0, so2],
            met={"TC0": np.array([295.0, 300.0]), ...},
            land_cover=[patches_cell0, patches_cell1],
            coefficients=olson_land_cover.coefficients,
        )
        # result → {"Hg0": array([...]), "SO2": array([...])}

    Parameters
    ----------
    species : list of mc.Species
        Each species must carry ``other_properties["henrys_law_constant"]``
        (M atm⁻¹) and ``other_properties["reactivity"]`` (0–1) as strings.
    box_height_m : float
        Height of the well-mixed box (m).
    met : dict
        Meteorological parameters; values may be scalars or 1-D arrays of
        length N.
    land_cover : list of LandCoverPatch or list of lists
        Flat list for scalar/single-cell mode; list of N patch lists for
        vectorized mode.
    coefficients : DrydepCoefficients
        Algorithm coefficients. Use ``olson_land_cover.coefficients`` for the
        standard Olson 2001 values.

    Returns
    -------
    single-species scalar mode
        ``(dvel_cms: float, k_s: float)``
    vectorized mode
        ``dict[str, np.ndarray]`` mapping species name → k array of length N.
    """
    _first_val = next(iter(met.values()))
    is_scalar = not hasattr(_first_val, "__len__")

    if len(species) == 1 and is_scalar and isinstance(land_cover[0], LandCoverPatch):
        sp = species[0]
        dvel_cms, k_s = _compute_single_cell(
            met, land_cover, sp, box_height_m, coefficients
        )
        return dvel_cms, k_s

    # Vectorized path
    if hasattr(_first_val, "__len__"):
        N = len(_first_val)
    else:
        N = len(land_cover)

    if isinstance(land_cover[0], LandCoverPatch):
        cell_land_cover = [land_cover] * N
    else:
        cell_land_cover = land_cover

    met_arrays = {}
    for key, val in met.items():
        if hasattr(val, "__len__"):
            met_arrays[key] = np.asarray(val, dtype=float)
        else:
            met_arrays[key] = np.full(N, float(val))

    result = {sp.name: np.zeros(N) for sp in species}

    for i in range(N):
        cell_met = {k: float(v[i]) for k, v in met_arrays.items()}
        patches = cell_land_cover[i]
        for sp in species:
            dvel_cms, k_s = _compute_single_cell(
                cell_met, patches, sp, box_height_m, coefficients
            )
            result[sp.name][i] = k_s

    return result


def _compute_single_cell(
    m: dict,
    patches: List[LandCoverPatch],
    sp: mc.Species,
    box_height_m: float,
    coefficients: DrydepCoefficients,
) -> tuple:
    """Run METERO + DEPVEL for one grid cell and one species."""
    CZ1, LSNOW, OBK, _ = METERO(
        m["BXHEIGHT"], m["ALBD"], m["TC0"], m["USTAR"],
        m["AIRDEN"], m["HFLUX"], m["U10M"], m["V10M"],
    )

    (IREG, ILAND, IUSE, XLAI,
     idep_arr, iri_arr, irlu_arr, irac_arr,
     irgss_arr, irgso_arr, ircls_arr, irclo_arr) = _patches_to_depvel_arrays(patches)

    hstar = float(sp.other_properties["henrys_law_constant"])
    f0 = float(sp.other_properties["reactivity"])

    dvel_cms = DEPVEL(
        coefficients.drycoeff, coefficients.iolson,
        idep_arr, iri_arr, irlu_arr, irac_arr,
        irgss_arr, irgso_arr, ircls_arr, irclo_arr,
        IREG, ILAND, IUSE, m["TC0"],
        XLAI, LSNOW, m["RADIAT"], m["CFRAC"], m["SUNCOS_MID"],
        m["PRESSU"], m["USTAR"], m["AZO"], CZ1, OBK,
        sp.molecular_weight_kg_mol, f0, hstar,
    )

    k_s = dvel_cms * 0.01 / box_height_m
    return dvel_cms, k_s
