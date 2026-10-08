"""
Aerosol Wet Deposition: CAM wetdepa_v2 (bulk / non-cloudborne branch)
=====================================================================
Target Parameterization for QUACS / CheMPAS-A

Python port of the CAM (CESM) aerosol wet-scavenging routines in
    src/chemistry/aerosol/wetdep.F90
        subroutine clddiag     -> cloud_volume_diag()
        subroutine wetdepa_v2  -> wetdepa_v2()    (is_strat_cloudborne = .false.,
                                                   no f_act_conv / qqcw)
    https://github.com/ESCOMP/CAM

References:
    Rasch, P. J., Feichter, J., Law, K., et al. (2000). A comparison of
    scavenging and deposition processes in global models: results from the
    WCRP Cambridge Workshop of 1995. Tellus B, 52, 1025-1056.
    https://doi.org/10.3402/tellusb.v52i4.17085

    Barth, M. C., Rasch, P. J., Kiehl, J. T., Benkovitz, C. M., and
    Schwartz, S. E. (2000). Sulfur chemistry in the National Center for
    Atmospheric Research Community Climate Model: Description, evaluation,
    features, and sensitivity to aqueous chemistry.
    J. Geophys. Res., 105, 1387-1415. https://doi.org/10.1029/1999JD900773

    Dana, M. T. and Hales, J. M. (1976). Statistical aspects of the washout
    of polydisperse aerosols. Atmos. Environ., 10, 45-50.
    https://doi.org/10.1016/0004-6981(76)90258-1

QUACS Contract
--------------
Standalone Python module. No CAM, no config files, no external dependencies
beyond NumPy. The tests are in tests/test_wetdep_cam.py. The examples in
examples/ run the scheme in a MUSICA box model and a MUSICA column model.

Array layout follows CAM: [ncol, nlev], level index k = 0 is the MODEL TOP
and k = nlev-1 is the level next to the surface.  All column sweeps go
top -> bottom, as in the Fortran.

Physics summary (wetdepa_v2, non-cloudborne branch), per level k:

  Precipitation and evaporation (fluxes accumulate top -> bottom):
      pdog      = pdel / g                                   [kg m-2]
      precabs   = sum_{above} (precs - evaps)  * pdog        [kg m-2 s-1]  stratiform
      precabc   = sum_{above} (cmfdqr - evapc) * pdog        [kg m-2 s-1]  convective
      fracev    = clip(evaps * pdog / precabs, 0, 1)         [-]  frac. of strat precip
                                                                  evaporating in layer k
      fracev_cu = clip(evapc * pdog / precabc, 0, 1)         [-]  same, convective

  In-cloud (nucleation) scavenging — fraction of cloud water converted
  to precipitation in one time step:
      stratiform:  fracp = clip(precs*dt / (cwat + precs*dt), 0, 1)
                   st_ic = sol_facti  * clds * fracp * q / dt     [kg kg-1 s-1]
      convective:  fracp = clip(cmfdqr*dt / (cldc*conicw + (cmfdqr+dlf)*dt), 0, 1)
                   cv_ic = sol_factic * cldc * fracp * q / dt     [kg kg-1 s-1]
      (clds = cldt - cldc)

  Below-cloud (impaction) scavenging, Dana & Hales (1976) coefficient:
      odds_s = clip(precabs / max(cldvst,1e-5) * scavcoef * dt, 0, 1)   [-]
      st_bc  = sol_factb * cldvst * odds_s * q / dt                     [kg kg-1 s-1]
      odds_c = clip(precabc / max(cldvcu,1e-5) * scavcoef * dt, 0, 1)
      cv_bc  = sol_factb * cldvcu * odds_c * q / dt
      Units of odds: [kg m-2 s-1] * [mm-1] * [s]; 1 mm of water = 1 kg m-2,
      so scavcoef [mm-1] = [m2 kg-1] and odds is dimensionless.
      (precabs / cldvst is the precipitation rate *inside* the raining area.)

  Totals and limiter (never remove more than q in one step):
      srcs = st_ic + st_bc ,  srcc = cv_ic + cv_bc
      rat  = q / (dt * (srcs + srcc));  if rat < 1: srcs, srcc *= rat
      srct = (srcs + srcc) * omsm                 omsm = 1 - 2*eps

  Tendency including resuspension of mass scavenged above:
      scavt = -srct + (fracev*scavab + fracev_cu*scavabc) * g / pdel    [kg kg-1 s-1]
      scavab  <- scavab *(1-fracev)    + srcs * pdog                   [kg m-2 s-1]
      scavabc <- scavabc*(1-fracev_cu) + srcc * pdog
  Surface wet deposition flux = scavab + scavabc below the lowest level
                                = -sum_k scavt_k * pdog_k               [kg m-2 s-1]

  Because every removal term is proportional to q, the removal is a
  first-order loss with rate coefficients (independent of q)
      k_strat = srcs/q ,  k_conv = srcc/q      [s-1]   (after the limiter)
  which is what the MUSICA box/column models hand to MICM as LOSS rates.

Cloud "volume" (precipitating area) diagnostic, clddiag (top -> bottom):
      lprec_st  = pdog * (precs - evaps);   lprecp = max(lprec_st, 1e-30)
      cldvst(k) = clip(cldv1_st/sumpppr_st, 0, 1) * (sumppr_st/sumpppr_st)
                  evaluated with sums over levels ABOVE k, then
      cldv1_st += clds*lprecp;  sumppr_st += lprec_st;  sumpppr_st += lprecp
      (cldvcu analogous with cmfdqr, evapc, cldc)
  i.e. the precipitation-weighted cloud fraction of the layers above.

Inputs (per time step, from the host model):
    q        : [ncol, nlev]  tracer mass mixing ratio                 [kg kg-1]
    pdel     : [ncol, nlev]  layer pressure thickness                 [Pa]
    cldt     : [ncol, nlev]  total cloud fraction                     [-]
    cldc     : [ncol, nlev]  convective cloud fraction                [-]
    cwat     : [ncol, nlev]  stratiform cloud water                   [kg kg-1]
    precs    : [ncol, nlev]  stratiform precip production rate        [kg kg-1 s-1]
    evaps    : [ncol, nlev]  stratiform precip evaporation rate       [kg kg-1 s-1]
    conicw   : [ncol, nlev]  in-cloud convective water                [kg kg-1]
    cmfdqr   : [ncol, nlev]  convective precip production rate        [kg kg-1 s-1]
    evapc    : [ncol, nlev]  convective precip evaporation rate       [kg kg-1 s-1]
    dlf      : [ncol, nlev]  detrained convective water (optional)    [kg kg-1 s-1]

Configurable (tracer-specific in CAM; set by the calling aerosol module):
    sol_fact  : solubility factor (fraction of tracer available to
                scavenging).  sol_facti (stratiform in-cloud),
                sol_factic (convective in-cloud) and sol_factb
                (below-cloud) default to sol_fact, as in wetdepa_v2.
    scavcoef  : below-cloud (Dana & Hales) coefficient [mm-1].
                CAM comment: "0.1 if not MODAL_AERO".  For MAM it is
                size-dependent (modal_aero_bcscavcoef); not ported here.

Output:
    scavt     : [ncol, nlev]  tracer tendency                         [kg kg-1 s-1]
    sfc_flux  : [ncol]        surface wet deposition flux             [kg m-2 s-1]

Not ported (out of scope for a bulk/non-modal tracer):
    cloud-borne branch (is_strat_cloudborne, qqcw, f_act_conv),
    MAM size-dependent below-cloud coefficients, ice/snow-specific factors.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Constants (CAM / shr_const)
# ---------------------------------------------------------------------------
GRAVIT = 9.80616  # [m s-2]
OMSM = 1.0 - 2.0 * np.finfo(float).eps  # roundoff guard (CAM omsm)
SCAVCOEF_DEFAULT = 0.1  # [mm-1] Dana & Hales, non-modal CAM
CLDV_MIN = 1.0e-5  # [-]   floor on raining area in odds
PREC_MIN = 1.0e-12  # [kg m-2 s-1] / [kg kg-1] floors (CAM)
TINY = 1.0e-36  # CAM 1.e-36_r8


# ---------------------------------------------------------------------------
# clddiag: precipitating cloud "volume" (area) seen by below-cloud scavenging
# ---------------------------------------------------------------------------
def cloud_volume_diag(pdel, cldt, cldc, precs, evaps, cmfdqr, evapc):
    """
    Port of CAM clddiag (wetdep.F90): precipitation-weighted cloud fraction
    of the layers ABOVE each level, used as the raining area for below-cloud
    scavenging.

        lprec_st   = (pdel/g) * (precs - evaps)                [kg m-2 s-1]
        lprecp_st  = max(lprec_st, 1e-30)
        cldvst(k)  = clip(cldv1_st/sumpppr_st, 0, 1) * sumppr_st/sumpppr_st,
                     floored at 0, with sums over levels above k
        cldv1_st  += (cldt - cldc) * lprecp_st
        sumppr_st += lprec_st;   sumpppr_st += lprecp_st

    cldvcu is the same with (cmfdqr, evapc, cldc).
    The sums start at 0 (sumpppr at TINY to avoid 0/0), so cldv = 0 at the
    model top.

    Parameters
    ----------
    pdel, cldt, cldc, precs, evaps, cmfdqr, evapc : np.ndarray [ncol, nlev]
        see module docstring for units.

    Returns
    -------
    cldvst : np.ndarray [ncol, nlev]  stratiform raining area [-]
    cldvcu : np.ndarray [ncol, nlev]  convective raining area [-]
    """
    ncol, nlev = pdel.shape
    clds = cldt - cldc

    cldv1_st = np.zeros(ncol)
    sumppr_st = np.zeros(ncol)
    sumpppr_st = np.full(ncol, TINY)
    cldv1_cu = np.zeros(ncol)
    sumppr_cu = np.zeros(ncol)
    sumpppr_cu = np.full(ncol, TINY)
    cldvst = np.zeros((ncol, nlev))
    cldvcu = np.zeros((ncol, nlev))

    for k in range(nlev):
        pdog = pdel[:, k] / GRAVIT

        # --- stratiform -----------------------------------------------------
        cldvst[:, k] = np.maximum(
            np.minimum(1.0, cldv1_st / sumpppr_st) * (sumppr_st / sumpppr_st), 0.0
        )
        lprec_st = pdog * (precs[:, k] - evaps[:, k])
        lprecp_st = np.maximum(lprec_st, 1.0e-30)
        cldv1_st += clds[:, k] * lprecp_st
        sumppr_st += lprec_st
        sumpppr_st += lprecp_st

        # --- convective -----------------------------------------------------
        cldvcu[:, k] = np.maximum(
            np.minimum(1.0, cldv1_cu / sumpppr_cu) * (sumppr_cu / sumpppr_cu), 0.0
        )
        lprec_cu = pdog * (cmfdqr[:, k] - evapc[:, k])
        lprecp_cu = np.maximum(lprec_cu, 1.0e-30)
        cldv1_cu += cldc[:, k] * lprecp_cu
        sumppr_cu += lprec_cu
        sumpppr_cu += lprecp_cu

    return cldvst, cldvcu


# ---------------------------------------------------------------------------
# Single-layer first-order scavenging rates (the core of wetdepa_v2)
# ---------------------------------------------------------------------------
def layer_scavenging_rates(
    cldt,
    cldc,
    cwat,
    precs,
    conicw,
    cmfdqr,
    precabs,
    precabc,
    cldvst,
    cldvcu,
    deltat,
    sol_fact,
    scavcoef=SCAVCOEF_DEFAULT,
    sol_facti=None,
    sol_factic=None,
    dlf=None,
):
    """
    First-order wet-scavenging rate coefficients for ONE layer, i.e. the
    wetdepa_v2 removal terms divided by the tracer mixing ratio q.

        k_st_ic = sol_facti  * clds   * fracp_s / dt
        k_st_bc = sol_factb  * cldvst * odds_s  / dt
        k_cv_ic = sol_factic * cldc   * fracp_c / dt
        k_cv_bc = sol_factb  * cldvcu * odds_c  / dt
        k_strat = k_st_ic + k_st_bc,  k_conv = k_cv_ic + k_cv_bc     [s-1]

    then the wetdepa_v2 limiter (no more than q removed per step):
        if dt*(k_strat + k_conv) > 1: scale both by 1/(dt*(k_strat + k_conv))
    and the roundoff guard omsm is applied to both.

    Parameters
    ----------
    cldt, cldc  : total / convective cloud fraction                     [-]
    cwat        : stratiform cloud water                                [kg kg-1]
    precs       : stratiform precip production                          [kg kg-1 s-1]
    conicw      : in-cloud convective water                             [kg kg-1]
    cmfdqr      : convective precip production                          [kg kg-1 s-1]
    precabs     : stratiform precip flux from layers above              [kg m-2 s-1]
    precabc     : convective precip flux from layers above              [kg m-2 s-1]
    cldvst      : stratiform raining area (clddiag)                     [-]
    cldvcu      : convective raining area (clddiag)                     [-]
    deltat      : time step                                             [s]
    sol_fact    : solubility factor (default for all three factors)     [-]
    scavcoef    : below-cloud coefficient                               [mm-1]
    sol_facti   : stratiform in-cloud factor (default sol_fact)         [-]
    sol_factic  : convective in-cloud factor (default sol_facti)        [-]
    dlf         : detrained convective water (default 0)                [kg kg-1 s-1]

    Returns
    -------
    dict of np.ndarray (same shape as inputs), all [s-1] except 'fracp_*'/'odds_*' [-]:
        k_strat, k_conv, k_total     rates after limiter and omsm
        k_st_ic, k_st_bc, k_cv_ic, k_cv_bc   components before the limiter
        fin_strat, fin_conv          in-cloud share of each process [-]
    """
    cldt = np.asarray(cldt, dtype=float)
    shape = np.broadcast(
        cldt, cldc, cwat, precs, conicw, cmfdqr, precabs, precabc, cldvst, cldvcu
    ).shape
    if dlf is None:
        dlf = np.zeros(shape)
    sol_facti = sol_fact if sol_facti is None else sol_facti
    sol_factic = sol_facti if sol_factic is None else sol_factic
    sol_factb = sol_fact
    dt = float(deltat)
    clds = cldt - cldc

    # --- convective ---------------------------------------------------------
    fracp_c = cmfdqr * dt / np.maximum(PREC_MIN, cldc * conicw + (cmfdqr + dlf) * dt)
    fracp_c = np.clip(fracp_c, 0.0, 1.0)
    k_cv_ic = sol_factic * cldc * fracp_c / dt
    odds_c = np.clip(precabc / np.maximum(cldvcu, CLDV_MIN) * scavcoef * dt, 0.0, 1.0)
    k_cv_bc = sol_factb * cldvcu * odds_c / dt

    # --- stratiform ---------------------------------------------------------
    fracp_s = precs * dt / np.maximum(PREC_MIN, cwat + precs * dt)
    fracp_s = np.clip(fracp_s, 0.0, 1.0)
    k_st_ic = sol_facti * clds * fracp_s / dt
    odds_s = np.clip(precabs / np.maximum(cldvst, CLDV_MIN) * scavcoef * dt, 0.0, 1.0)
    k_st_bc = sol_factb * cldvst * odds_s / dt

    k_conv = k_cv_ic + k_cv_bc
    k_strat = k_st_ic + k_st_bc
    fin_conv = k_cv_ic / (k_conv + TINY)
    fin_strat = k_st_ic / (k_strat + TINY)

    # --- limiter: rat = q / (dt*(srcc+srcs)) with q factored out -------------
    rat = 1.0 / np.maximum(dt * (k_conv + k_strat), TINY)
    scale = np.where(rat < 1.0, rat, 1.0)
    k_conv = k_conv * scale * OMSM
    k_strat = k_strat * scale * OMSM

    return dict(
        k_strat=k_strat,
        k_conv=k_conv,
        k_total=k_strat + k_conv,
        k_st_ic=k_st_ic,
        k_st_bc=k_st_bc,
        k_cv_ic=k_cv_ic,
        k_cv_bc=k_cv_bc,
        fin_strat=fin_strat,
        fin_conv=fin_conv,
        fracp_s=fracp_s,
        fracp_c=fracp_c,
        odds_s=odds_s,
        odds_c=odds_c,
    )


# ---------------------------------------------------------------------------
# Column coefficients (tracer-independent): rates + evaporation fractions
# ---------------------------------------------------------------------------
def column_scavenging_coefficients(
    pdel,
    cldt,
    cldc,
    cwat,
    precs,
    evaps,
    conicw,
    cmfdqr,
    evapc,
    deltat,
    sol_fact,
    scavcoef=SCAVCOEF_DEFAULT,
    sol_facti=None,
    sol_factic=None,
    dlf=None,
    cldvst=None,
    cldvcu=None,
):
    """
    Sweep a column top -> bottom and return the tracer-independent pieces of
    wetdepa_v2: first-order rates (k_strat, k_conv) and the fractions of the
    precipitation flux from above that evaporates in each layer
    (fracev, fracev_cu), which control resuspension.

    Returns
    -------
    dict of np.ndarray [ncol, nlev]:
        k_strat, k_conv   [s-1]
        fracev, fracev_cu [-]
        pdog              [kg m-2]
        precabs, precabc  [kg m-2 s-1]  precipitation flux entering each layer from above
        cldvst, cldvcu    [-]
        plus the per-process components (k_st_ic, k_st_bc, k_cv_ic, k_cv_bc)
    """
    ncol, nlev = pdel.shape
    if dlf is None:
        dlf = np.zeros_like(pdel)
    if cldvst is None or cldvcu is None:
        cldvst, cldvcu = cloud_volume_diag(pdel, cldt, cldc, precs, evaps, cmfdqr, evapc)

    out = {
        name: np.zeros((ncol, nlev))
        for name in (
            "k_strat",
            "k_conv",
            "fracev",
            "fracev_cu",
            "pdog",
            "precabs",
            "precabc",
            "k_st_ic",
            "k_st_bc",
            "k_cv_ic",
            "k_cv_bc",
        )
    }
    out["cldvst"], out["cldvcu"] = cldvst, cldvcu

    precabs = np.zeros(ncol)  # [kg m-2 s-1] stratiform precip flux from above
    precabc = np.zeros(ncol)  # [kg m-2 s-1] convective precip flux from above

    for k in range(nlev):
        pdog = pdel[:, k] / GRAVIT  # [kg m-2]
        fracev = np.clip(evaps[:, k] * pdog / np.maximum(PREC_MIN, precabs), 0.0, 1.0)
        fracev_cu = np.clip(evapc[:, k] * pdog / np.maximum(PREC_MIN, precabc), 0.0, 1.0)

        r = layer_scavenging_rates(
            cldt[:, k],
            cldc[:, k],
            cwat[:, k],
            precs[:, k],
            conicw[:, k],
            cmfdqr[:, k],
            precabs,
            precabc,
            cldvst[:, k],
            cldvcu[:, k],
            deltat,
            sol_fact,
            scavcoef,
            sol_facti,
            sol_factic,
            dlf[:, k],
        )

        out["k_strat"][:, k] = r["k_strat"]
        out["k_conv"][:, k] = r["k_conv"]
        out["k_st_ic"][:, k] = r["k_st_ic"]
        out["k_st_bc"][:, k] = r["k_st_bc"]
        out["k_cv_ic"][:, k] = r["k_cv_ic"]
        out["k_cv_bc"][:, k] = r["k_cv_bc"]
        out["fracev"][:, k] = fracev
        out["fracev_cu"][:, k] = fracev_cu
        out["pdog"][:, k] = pdog
        out["precabs"][:, k] = precabs
        out["precabc"][:, k] = precabc

        # precipitation flux leaving the bottom of layer k (CAM order: after use)
        precabs = precabs + (precs[:, k] - evaps[:, k]) * pdog
        precabc = precabc + (cmfdqr[:, k] - evapc[:, k]) * pdog

    return out


# ---------------------------------------------------------------------------
# wetdepa_v2: tendency with resuspension (CAM explicit form)
# ---------------------------------------------------------------------------
def wetdepa_v2(
    q,
    pdel,
    cldt,
    cldc,
    cwat,
    precs,
    evaps,
    conicw,
    cmfdqr,
    evapc,
    deltat,
    sol_fact,
    scavcoef=SCAVCOEF_DEFAULT,
    sol_facti=None,
    sol_factic=None,
    dlf=None,
    cldvst=None,
    cldvcu=None,
):
    """
    CAM wetdepa_v2 (non-cloudborne branch): tracer tendency due to wet
    scavenging, including resuspension where precipitation evaporates.

        srcs  = k_strat * q,  srcc = k_conv * q                     [kg kg-1 s-1]
        scavt = -(srcs + srcc) + (fracev*scavab + fracev_cu*scavabc) * g/pdel
        scavab  <- scavab *(1-fracev)    + srcs * pdel/g            [kg m-2 s-1]
        scavabc <- scavabc*(1-fracev_cu) + srcc * pdel/g

    (k_* already include the limiter and omsm; in CAM omsm multiplies srct
    but not the srcs added to scavab — the difference is O(1e-16).)

    CAM applies this explicitly:  q_new = q + scavt * deltat.

    Parameters
    ----------
    q : np.ndarray [ncol, nlev]  tracer mass mixing ratio [kg kg-1]
    (others: see column_scavenging_coefficients)

    Returns
    -------
    dict:
        scavt     [ncol, nlev]  total tendency                      [kg kg-1 s-1]
        ic_scavt  [ncol, nlev]  in-cloud removal part (<= 0)        [kg kg-1 s-1]
        bc_scavt  [ncol, nlev]  below-cloud removal part (<= 0)     [kg kg-1 s-1]
        resusp    [ncol, nlev]  resuspension source (>= 0)          [kg kg-1 s-1]
        sfc_flux  [ncol]        surface wet deposition flux         [kg m-2 s-1]
        coef      dict from column_scavenging_coefficients
    """
    c = column_scavenging_coefficients(
        pdel,
        cldt,
        cldc,
        cwat,
        precs,
        evaps,
        conicw,
        cmfdqr,
        evapc,
        deltat,
        sol_fact,
        scavcoef,
        sol_facti,
        sol_factic,
        dlf,
        cldvst,
        cldvcu,
    )
    ncol, nlev = q.shape

    # in-cloud share of each process, recomputed from the pre-limiter components
    fin_s = c["k_st_ic"] / (c["k_st_ic"] + c["k_st_bc"] + TINY)
    fin_c = c["k_cv_ic"] / (c["k_cv_ic"] + c["k_cv_bc"] + TINY)

    scavt = np.zeros((ncol, nlev))
    ic_scavt = np.zeros((ncol, nlev))
    bc_scavt = np.zeros((ncol, nlev))
    resusp = np.zeros((ncol, nlev))
    scavab = np.zeros(ncol)  # [kg m-2 s-1] tracer flux in strat precip from above
    scavabc = np.zeros(ncol)  # [kg m-2 s-1] tracer flux in conv precip from above

    for k in range(nlev):
        pdog = c["pdog"][:, k]
        srcs = c["k_strat"][:, k] * q[:, k]
        srcc = c["k_conv"][:, k] * q[:, k]
        res = (c["fracev"][:, k] * scavab + c["fracev_cu"][:, k] * scavabc) / pdog

        scavt[:, k] = -(srcs + srcc) + res
        ic_scavt[:, k] = -(srcs * fin_s[:, k] + srcc * fin_c[:, k])
        bc_scavt[:, k] = -(srcs * (1.0 - fin_s[:, k]) + srcc * (1.0 - fin_c[:, k]))
        resusp[:, k] = res

        scavab = scavab * (1.0 - c["fracev"][:, k]) + srcs * pdog
        scavabc = scavabc * (1.0 - c["fracev_cu"][:, k]) + srcc * pdog

    return dict(
        scavt=scavt,
        ic_scavt=ic_scavt,
        bc_scavt=bc_scavt,
        resusp=resusp,
        sfc_flux=scavab + scavabc,
        coef=c,
    )


# ---------------------------------------------------------------------------
# Exact-in-time column step (used to drive / verify MICM)
# ---------------------------------------------------------------------------
def column_step_exponential(q, coef, deltat):
    """
    Advance one time step with the same top -> bottom resuspension recursion
    as wetdepa_v2, but integrating each layer EXACTLY over the step:

        dq/dt = -k q + s,   k = k_strat + k_conv,   s = resuspension source
        q_new = s/k + (q - s/k) * exp(-k dt)           (k > 0)
        q_new = q + s dt                               (k = 0)

    The resuspension source s of layer k is built from the mass actually
    removed (step-mean) in the layers above, so column mass is conserved
    exactly:  burden_old - burden_new = sfc_flux * dt.

    For k*dt << 1 this equals CAM's explicit update q + scavt*dt.

    Parameters
    ----------
    q      : np.ndarray [ncol, nlev]  [kg kg-1]
    coef   : dict from column_scavenging_coefficients
    deltat : float [s]

    Returns
    -------
    q_new    : [ncol, nlev]  [kg kg-1]
    s        : [ncol, nlev]  resuspension source, constant over the step [kg kg-1 s-1]
    sfc_flux : [ncol]        step-mean surface wet deposition flux     [kg m-2 s-1]
    """
    ncol, nlev = q.shape
    dt = float(deltat)
    q_new = np.zeros_like(q)
    s_all = np.zeros_like(q)
    scavab = np.zeros(ncol)
    scavabc = np.zeros(ncol)

    for k in range(nlev):
        pdog = coef["pdog"][:, k]
        ks, kc = coef["k_strat"][:, k], coef["k_conv"][:, k]
        kt = ks + kc
        s = (coef["fracev"][:, k] * scavab + coef["fracev_cu"][:, k] * scavabc) / pdog

        e = np.exp(-kt * dt)
        with np.errstate(divide="ignore", invalid="ignore"):
            qn = np.where(kt > 0.0, s / kt + (q[:, k] - s / kt) * e, q[:, k] + s * dt)
        # step-mean removal rate R = (1/dt) * int k q dt = (q - q_new)/dt + s
        R = (q[:, k] - qn) / dt + s
        fs = np.where(kt > 0.0, ks / np.where(kt > 0.0, kt, 1.0), 0.0)

        scavab = scavab * (1.0 - coef["fracev"][:, k]) + R * fs * pdog
        scavabc = scavabc * (1.0 - coef["fracev_cu"][:, k]) + R * (1.0 - fs) * pdog
        q_new[:, k] = qn
        s_all[:, k] = s

    return q_new, s_all, scavab + scavabc


def column_burden(q, pdel):
    """Column-integrated tracer mass  sum_k q * pdel / g   [kg m-2]."""
    return (q * pdel / GRAVIT).sum(axis=-1)
