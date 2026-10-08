"""Print the CAM wetdepa_v2 tendency for one column with one stratiform cloud layer."""

import numpy as np

from quacs.wetdep.aerosol_cam import GRAVIT, wetdepa_v2


def _single_cloud_column(
    nlev=6, cloud_k=1, cldfrac=0.4, precs_val=1.0e-7, evap_k=None, evap_frac=0.0
):
    """Helper: one stratiform cloud layer at cloud_k, optional evaporation layer."""
    shp = (1, nlev)
    pdel = np.full(shp, 5000.0)
    cldt = np.zeros(shp)
    cldt[0, cloud_k] = cldfrac
    cldc = np.zeros(shp)
    cwat = np.zeros(shp)
    cwat[0, cloud_k] = 2.0e-4
    precs = np.zeros(shp)
    precs[0, cloud_k] = precs_val
    evaps = np.zeros(shp)
    if evap_k is not None:
        flux = precs_val * 5000.0 / GRAVIT
        evaps[0, evap_k] = evap_frac * flux * GRAVIT / 5000.0
    zeros = np.zeros(shp)
    return dict(
        pdel=pdel,
        cldt=cldt,
        cldc=cldc,
        cwat=cwat,
        precs=precs,
        evaps=evaps,
        conicw=zeros,
        cmfdqr=zeros,
        evapc=zeros,
    )


col = _single_cloud_column(nlev=6, cloud_k=1, evap_k=4, evap_frac=0.4)
q = np.full((1, 6), 3.0e-8)
r = wetdepa_v2(q=q, deltat=1800.0, sol_fact=0.3, **col)
print("CAM wetdepa_v2 demo (one stratiform cloud layer at k=1, evaporation at k=4)")
print(
    f"{'k':>3} {'k_strat [s-1]':>15} {'cldvst':>8} {'precabs':>12} "
    f"{'scavt [kg/kg/s]':>17} {'resusp':>12}"
)
for k in range(6):
    print(
        f"{k:>3} {r['coef']['k_strat'][0, k]:>15.4e} {r['coef']['cldvst'][0, k]:>8.3f} "
        f"{r['coef']['precabs'][0, k]:>12.4e} {r['scavt'][0, k]:>17.4e} "
        f"{r['resusp'][0, k]:>12.4e}"
    )
print(f"surface wet deposition flux = {r['sfc_flux'][0]:.4e} kg m-2 s-1")
