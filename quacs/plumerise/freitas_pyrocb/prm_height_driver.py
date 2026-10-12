"""Step 2 driver: resolve baseline and pyroCb injection heights."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

INJECTH_BASE = 3000.0
INJECTH_TOP = 6000.0
PYROCB_HALF_WIDTH = 2000.0
FALLBACK_TROPOPAUSE_HEIGHT = 12000.0


@dataclass(frozen=True)
class TropopauseResult:
    tropopause_height_m: Optional[float]
    source: str
    status_flag: str
    diagnostics: dict = field(default_factory=dict)


def compute_idealized_prm_heights(
    z,
    p,
    t,
    u,
    v,
    qv,
    vegetation_class,
    fire_size_mean,
    fire_size_std,
):
    _ = z, p, t, u, v, qv, vegetation_class, fire_size_mean, fire_size_std
    return INJECTH_BASE, INJECTH_TOP


def compute_wmo_lapse_rate_tropopause(
    z,
    t,
):
    if len(z) != len(t) or len(z) < 3:
        return TropopauseResult(None, "missing", "invalid_profile")
    if np.any(np.diff(z) <= 0.0):
        return TropopauseResult(None, "missing", "invalid_profile")

    for idx in range(len(z) - 1):
        dz_next = z[idx + 1] - z[idx]
        local_lapse = -(t[idx + 1] - t[idx]) / dz_next * 1000.0
        if local_lapse > 2.0:
            continue
        average_lapses = []
        for jdx in range(idx + 1, len(z)):
            if z[jdx] - z[idx] > 2000.0:
                break
            average_lapses.append(-(t[jdx] - t[idx]) / (z[jdx] - z[idx]) * 1000.0)
        if average_lapses and max(average_lapses) <= 2.0:
            return TropopauseResult(
                tropopause_height_m=float(z[idx]),
                source="wmo",
                status_flag="ok",
                diagnostics={
                    "local_lapse_K_km": float(local_lapse),
                    "max_average_lapse_K_km": float(max(average_lapses)),
                },
            )

    return TropopauseResult(
        tropopause_height_m=FALLBACK_TROPOPAUSE_HEIGHT,
        source="fallback",
        status_flag="tropopause_fallback_used",
        diagnostics={"fallback_tropopause_height": FALLBACK_TROPOPAUSE_HEIGHT},
    )


def compute_pyrocb_height_range(
    pyrocb_flag,
    tropopause,
):
    if not pyrocb_flag or tropopause.tropopause_height_m is None:
        return None, None
    center = float(tropopause.tropopause_height_m)
    return (
        center - PYROCB_HALF_WIDTH,
        center + PYROCB_HALF_WIDTH,
    )


def prm_height_driver(
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
):
    injectH_base, injectH_top = compute_idealized_prm_heights(
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
    tropopause = compute_wmo_lapse_rate_tropopause(z, t)
    injectH_pyroCb_base, injectH_pyroCb_top = compute_pyrocb_height_range(
        pyrocb_flag,
        tropopause,
    )
    return (
        injectH_base,
        injectH_top,
        injectH_pyroCb_base,
        injectH_pyroCb_top,
    )
