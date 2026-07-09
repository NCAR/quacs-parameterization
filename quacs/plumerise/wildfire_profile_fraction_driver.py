"""Outer driver for the wildfire vertical profile fraction workflow."""

from __future__ import annotations

import numpy as np

from quacs.plumerise.layer_fraction_driver import layer_fraction_driver
from quacs.plumerise.prm_height_driver import prm_height_driver
from quacs.plumerise.pyrocb_flag_driver import pyrocb_flag_driver


def compute_wildfire_profile_fraction_driver(
    z,
    p,
    t,
    u,
    v,
    qv,
    vegetation_class,
    fire_size_mean,
    fire_size_std=0.0,
):
    """Compute the public `layer_fraction[k]` output from direct variables.

    Variables:
    - `z` [m AGL]: vertical coordinate for the profile. Step 3 derives layer
      bounds from `z` by appending one top bound using the final spacing.
    - `p` [hPa]: pressure profile on `z`; required by the pyrometeopy-style
      PFT calculation.
    - `t` [K]: temperature profile on `z`.
    - `u`, `v` [m s^-1]: horizontal wind profiles on `z`.
    - `qv` [kg kg^-1]: water-vapor specific humidity profile on `z`.
    - `vegetation_class` [unitless]: vegetation key used by the internal split
      and heat-flux tables.
    - `fire_size_mean` [m^2]: mean fire area / fire-size proxy for the profile.
    - `fire_size_std` [m^2]: optional fire-size uncertainty for the high-bound
      qplume estimate.

    Output:
    - `layer_fraction`: normalized vertical profile fraction for the input
      profile.
    """

    z = np.asarray(z, dtype=float)
    p = np.asarray(p, dtype=float)
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    qv = np.asarray(qv, dtype=float)
    if len(z) < 2:
        raise ValueError("z must contain at least two levels")

    # Step 1: diagnose whether this profile should use the pyroCb branch.
    pyrocb_flag = pyrocb_flag_driver(
        z,
        p,
        t,
        u,
        v,
        qv,
        vegetation_class,
        fire_size_mean,
        fire_size_std,
    )

    # Step 2: resolve baseline PRM heights and any pyroCb-adjusted height range.
    (
        injectH_base_m,
        injectH_top_m,
        injectH_pyroCb_base_m,
        injectH_pyroCb_top_m,
    ) = prm_height_driver(
        z,
        p,
        t,
        u,
        v,
        qv,
        vegetation_class,
        fire_size_mean,
        fire_size_std,
        pyrocb_flag,
    )

    # Step 3: convert direct weights and target height ranges to layer fractions.
    layer_fraction = layer_fraction_driver(
        z,
        vegetation_class,
        injectH_base_m,
        injectH_top_m,
        injectH_pyroCb_base_m,
        injectH_pyroCb_top_m,
        pyrocb_flag,
    )
    return layer_fraction
