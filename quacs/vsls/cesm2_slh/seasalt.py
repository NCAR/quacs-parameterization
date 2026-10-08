"""
Sea-salt bromine recycling rates following CESM2-SLH (CAM ``mo_usrrxt.F90``).

Computes pseudo-first-order rate constants for the heterogeneous uptake of
reactive bromine on sea-salt aerosol, which releases Br2 and BrCl to the gas
phase (sea-salt dehalogenation):

=========  ===========  ======================================  ======
CAM tag    Reactant     Products (``chem_mech.in``, SLH)        gamma
=========  ===========  ======================================  ======
het_ss_0   BrONO2       0.65 Br2 + 0.35 BrCl                    0.010
het_ss_1   BrNO2        0.65 Br2 + 0.35 BrCl                    0.005
het_ss_2   HOBr         0.65 Br2 + 0.35 BrCl                    0.0125
=========  ===========  ======================================  ======

Rate constant (free-molecular uptake, as implemented in CAM)::

    k = S * 0.25 * gamma * v_mean * A_ss * DF * mask      (1/s)

where ``v_mean = sqrt(8 R T / (pi M))`` (m/s), ``A_ss`` is the sea-salt
surface area density (m2/m3), ``DF`` is the bromine depletion factor
(Yang et al., 2005, as cited in CAM), ``S`` is the ``SSAdehal_ScalingFactor``
namelist value (CAM default 1.0), and ``mask`` reproduces CAM's
geographic/vertical masking.

Source (checked against ESCOMP/CAM ``cam_development`` branch):
    src/chemistry/mozart/mo_usrrxt.F90                 gamma, DF, masks, rate form
    src/chemistry/mozart/mo_slh_routines.F90           SSAdehal_ScalingFactor = 1.0
    src/chemistry/pp_trop_strat_mam4_slh/chem_mech.in  products and yields

Notes
-----
* Gas-phase bromine is intentionally not conserved: each reactant Br atom
  yields 0.65*2 + 0.35 = 1.65 gas-phase Br atoms. The extra 0.65 Br comes from
  sea-salt bromide, which is not tracked as a gas-phase species.
* CAM hard-codes the speed prefactor (e.g. 1.47e3 cm/s/K^0.5 for HOBr);
  here it is computed from molar mass. The two agree to <1% (see tests).
* CAM's free-molecular form omits gas-phase diffusion limitation. This is a
  CAM-consistent simplification, not a general uptake formula.
"""

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

R_GAS = 8.314462618  # J/mol/K

# Depletion-factor bounds south of 30 S (CAM: dfmax, dfmin)
DF_MAX = 0.9
DF_MIN = 0.3
DF_NORTH = 0.5  # constant DF north of 30 S

# Masking thresholds (CAM)
MIN_PRESSURE_PA = 300.0e2  # no sea-salt recycling above the 300 hPa level
LAND_MASK_LAT_DEG = -60.0  # pure-land points north of 60 S get zero sea-salt SAD


@dataclass(frozen=True)
class SeaSaltReaction:
    """One sea-salt dehalogenation reaction.

    Attributes
    ----------
    name : str
        Reaction identifier (CAM tag).
    reactant : str
        Gas-phase reactant species name.
    molar_mass_kg_mol : float
        Molar mass of the reactant (kg/mol), used for the mean speed.
    gamma : float
        Uptake coefficient (dimensionless).
    products : tuple of (str, float)
        Gas-phase products and molar yields.
    """

    name: str
    reactant: str
    molar_mass_kg_mol: float
    gamma: float
    products: Tuple[Tuple[str, float], ...]


_BR_YIELDS = (("Br2", 0.65), ("BrCl", 0.35))

SEA_SALT_REACTIONS: Tuple[SeaSaltReaction, ...] = (
    SeaSaltReaction("het_ss_0", "BrONO2", 141.91e-3, 0.010, _BR_YIELDS),
    SeaSaltReaction("het_ss_1", "BrNO2", 125.91e-3, 0.005, _BR_YIELDS),
    SeaSaltReaction("het_ss_2", "HOBr", 96.91e-3, 0.0125, _BR_YIELDS),
)


def mean_molecular_speed(temperature_k, molar_mass_kg_mol: float):
    """Mean molecular speed sqrt(8 R T / (pi M)) in m/s."""
    t = np.asarray(temperature_k, dtype=float)
    return np.sqrt(8.0 * R_GAS * t / (np.pi * molar_mass_kg_mol))


def bromine_depletion_factor(latitude_deg, calday):
    """Bromine depletion factor DF (dimensionless), CAM form.

    DF = 0.5 north of 30 S. South of 30 S::

        DF = DF_MAX + 0.5 (DF_MIN - DF_MAX) (sin((calday / 182.5 - 0.5) pi) + 1)

    Parameters
    ----------
    latitude_deg : float or array
        Latitude (degrees north).
    calday : float or array
        Calendar day of year, as used by CAM.
    """
    lat = np.asarray(latitude_deg, dtype=float)
    day = np.asarray(calday, dtype=float)
    south = DF_MAX + 0.5 * (DF_MIN - DF_MAX) * (np.sin((day / (365.0 / 2.0) - 0.5) * np.pi) + 1.0)
    return np.where(lat > -30.0, DF_NORTH, south)


def seasalt_mask(latitude_deg, ocean_ice_fraction, pressure_pa):
    """Return 1 where CAM applies sea-salt recycling, else 0.

    Recycling is off (a) over pure-land points north of 60 S, where CAM zeroes
    the sea-salt SAD, and (b) at pressures below 300 hPa.
    """
    lat = np.asarray(latitude_deg, dtype=float)
    ocn = np.asarray(ocean_ice_fraction, dtype=float)
    p = np.asarray(pressure_pa, dtype=float)
    land_off = (ocn <= 0.0) & (lat > LAND_MASK_LAT_DEG)
    upper_off = p < MIN_PRESSURE_PA
    return np.where(land_off | upper_off, 0.0, 1.0)


def seasalt_rate_constants(
    temperature_k,
    pressure_pa,
    latitude_deg,
    calday,
    sad_seasalt_m2_m3,
    ocean_ice_fraction,
    scaling_factor: float = 1.0,
) -> Dict[str, np.ndarray]:
    """Pseudo-first-order rate constants (1/s) for each sea-salt reaction.

    All inputs may be scalars or broadcast-compatible arrays.

    Parameters
    ----------
    temperature_k : float or array
        Air temperature (K).
    pressure_pa : float or array
        Mid-layer pressure (Pa).
    latitude_deg : float or array
        Latitude (degrees north).
    calday : float or array
        Calendar day of year (CAM convention).
    sad_seasalt_m2_m3 : float or array
        Sea-salt surface area density (m2/m3). CAM carries SAD in cm2/cm3;
        1 cm2/cm3 = 100 m2/m3.
    ocean_ice_fraction : float or array
        Ocean plus sea-ice fraction of the grid cell (0-1).
    scaling_factor : float
        CAM ``SSAdehal_ScalingFactor`` (default 1.0).

    Returns
    -------
    dict
        Reaction name -> rate constant array (1/s).
    """
    sad = np.asarray(sad_seasalt_m2_m3, dtype=float)
    if np.any(sad < 0.0):
        raise ValueError("sad_seasalt_m2_m3 must be non-negative")
    df = bromine_depletion_factor(latitude_deg, calday)
    mask = seasalt_mask(latitude_deg, ocean_ice_fraction, pressure_pa)
    rates = {}
    for rxn in SEA_SALT_REACTIONS:
        v = mean_molecular_speed(temperature_k, rxn.molar_mass_kg_mol)
        rates[rxn.name] = scaling_factor * 0.25 * rxn.gamma * v * sad * df * mask
    return rates
