"""MEGAN 3.0 biogenic VOC emissions for one grid cell.

>>> from quacs.biogenic.megan import MeganSettings, compute_emissions
"""

from .biogenic_emission_megan_v3 import MeganSettings
from .biogenic_emission_megan_v3 import biogenic_emission_megan_v3 as compute_emissions

__all__ = ["MeganSettings", "compute_emissions"]
