"""
GOCART dust emission scheme (Ginoux et al., 2001), as in MPAS-GOCART2G.

>>> from quacs.dust.ginoux import dust_emission, threshold_velocity, SP_DEFAULT
>>> from quacs.dust.ginoux import run_box_model
"""

from .box_model import build_mechanism, emission_rates, run_box_model
from .ginoux import (
    C_DEFAULT,
    GRAV,
    LAND,
    RADIUS_DEFAULT,
    RADIUS_LOWER,
    RADIUS_UPPER,
    RHOP_DEFAULT,
    SP_DEFAULT,
    dust_emission,
    initialize,
    threshold_velocity,
)

__all__ = [
    "C_DEFAULT",
    "GRAV",
    "LAND",
    "RADIUS_DEFAULT",
    "RADIUS_LOWER",
    "RADIUS_UPPER",
    "RHOP_DEFAULT",
    "SP_DEFAULT",
    "dust_emission",
    "initialize",
    "threshold_velocity",
    "build_mechanism",
    "emission_rates",
    "run_box_model",
]
