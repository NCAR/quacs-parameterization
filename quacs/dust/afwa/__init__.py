"""AFWA dust emission. See README.md for inputs and outputs."""

from .afwa import (
    AfwaSettings,
    afwa_dust_driver,
    dust_emission_afwa,
    initialize,
    saltation_area_fraction,
    threshold_friction_velocity,
    tunable_names,
)

__all__ = [
    "AfwaSettings",
    "afwa_dust_driver",
    "dust_emission_afwa",
    "initialize",
    "saltation_area_fraction",
    "threshold_friction_velocity",
    "tunable_names",
]
