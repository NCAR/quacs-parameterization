"""pyroCb-aware plume rise: the vertical profile fraction of wildfire emissions in one column.

The public driver ``compute_wildfire_profile_fraction_driver`` calls 3 steps:

1. ``pyrocb_flag_driver``: diagnose the pyroCb branch.
2. ``prm_height_driver``: resolve the injection height ranges.
3. ``layer_fraction_driver``: convert the height ranges to ``layer_fraction[k]``.
"""

from quacs.plumerise.freitas_pyrocb.layer_fraction_driver import layer_fraction_driver
from quacs.plumerise.freitas_pyrocb.prm_height_driver import prm_height_driver
from quacs.plumerise.freitas_pyrocb.pyrocb_flag_driver import pyrocb_flag_driver
from quacs.plumerise.freitas_pyrocb.wildfire_profile_fraction_driver import (
    compute_wildfire_profile_fraction_driver,
)

__all__ = [
    "compute_wildfire_profile_fraction_driver",
    "layer_fraction_driver",
    "prm_height_driver",
    "pyrocb_flag_driver",
]
