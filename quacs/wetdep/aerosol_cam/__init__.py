"""CAM aerosol wet deposition (wetdepa_v2, bulk, non-cloudborne branch).

>>> from quacs.wetdep.aerosol_cam import wetdepa_v2
>>> result = wetdepa_v2(q, pdel, cldt, cldc, cwat, precs, evaps, conicw, cmfdqr, evapc,
...                     deltat=1800.0, sol_fact=0.3)
>>> result["scavt"], result["sfc_flux"]
"""

from .wetdep_cam import (
    GRAVIT,
    SCAVCOEF_DEFAULT,
    cloud_volume_diag,
    column_burden,
    column_scavenging_coefficients,
    column_step_exponential,
    layer_scavenging_rates,
    wetdepa_v2,
)

__all__ = [
    "GRAVIT",
    "SCAVCOEF_DEFAULT",
    "cloud_volume_diag",
    "column_burden",
    "column_scavenging_coefficients",
    "column_step_exponential",
    "layer_scavenging_rates",
    "wetdepa_v2",
]
