"""
Dust Emission Scheme: DEAD / Zender et al. (2003a) — CLM4 implementation
=========================================================================
Target Parameterization for QUACS / CheMPAS-A

References:
    Zender, C. S., Bian, H. S., and Newman, D. (2003a). Mineral Dust
    Entrainment and Deposition (DEAD) model: Description and 1990s dust
    climatology. J. Geophys. Res.-Atmos., 108, 4416.
    https://doi.org/10.1029/2002JD002775

    Marticorena, B. and Bergametti, G. (1995). Modeling the atmospheric
    dust cycle: 1. Design of a soil-derived dust emission scheme.
    J. Geophys. Res., 100, 16415-16430.
    https://doi.org/10.1029/95JD00690

    Fécan, F., Marticorena, B., and Bergametti, G. (1999). Parametrization
    of the increase of the aeolian erosion threshold wind friction velocity
    due to soil moisture for arid and semi-arid areas. Ann. Geophys., 17,
    149-157. https://doi.org/10.1007/s00585-999-0149-7

    Kok, J. F., et al. (2014). An improved dust emission model – Part 2:
    Evaluation in the Community Earth System Model, with implications for
    the use of dust source functions. Atmos. Chem. Phys., 14, 13043-13061.
    https://doi.org/10.5194/acp-14-13043-2014

    Fortran reference: CLM / CTSM
      src/biogeochem/DustEmisZender2003.F90       (emission flux)
      src/biogeophys/SoilStateInitTimeConstMod.F90 (gwc_thr, clay cap)
    https://github.com/ESCOMP/CTSM

QUACS Contract
--------------
Standalone Python module. No MPAS-A, no config files, no external
dependencies beyond NumPy and SciPy. Import it as ``quacs.dust.zender``.
The unit tests are in ``tests/test_zender.py``. Run them with
``pytest quacs/dust/zender``. The examples are in ``examples/``.

Pre-API note:
    Once the QUACS scripting API is finalized, the expected signature is:
        @scheme(category='dust_emission')
        def dust_emission_dead(state: ModelState, config: Config) -> SurfaceFlux:
    Until then, plain NumPy arrays are used. This script is the scientific
    specification that drives that API design.

Physics summary (Zender et al. 2003a / Marticorena & Bergametti 1995):

  Step 1 — Mobilizable land fraction (fbare):
      fbare = max(0, 1 - VAI/vai_mbl_thr) * (1 - frac_sno)
      VAI   = tlai + tsai  (vegetation area index)
      Emission occurs only over soil/crop land types with fbare > 0.

  Step 2 — Soil moisture threshold factor (Fécan et al. 1999, Zender 2003a tuning):
      gwc_sfc = h2osoi_vol * rho_water / bd   [kg kg-1] gravimetric water content
      bd      = (1 - watsat) * 2700           [kg m-3]  bulk density of dry soil
      gwc_thr = 0.17 + 0.14*fclay             [kg kg-1] threshold gravimetric w.c.
      frc_thr_wet_fct = sqrt(1 + 1.21*(100*(gwc_sfc - gwc_thr))^0.68)  if gwc > gwc_thr
                      = 1.0                                              otherwise

      gwc_thr is Fécan's w' = 0.17*fclay + 0.14*fclay^2 multiplied by the
      Zender (2003a) tuning factor a = 1/fclay, i.e. 0.17 + 0.14*fclay.
      This is exactly CTSM's ThresholdSoilMoistZender2003
      (0.17 + 0.14*clay[%]*0.01).  The un-tuned Fécan form (a = 1) is what
      K14 uses, not DEAD/CLM.

  Step 3 — Threshold friction velocity for saltation (Iversen & White 1982
            via Marticorena & Bergametti 1995):
      tmp1 = precomputed constant for optimal saltation diameter (75 µm, 2650 kg/m³)
             [kg^0.5 m^-0.5 s^-1]
      wnd_frc_thr_slt = tmp1 / sqrt(forc_rho) * frc_thr_wet_fct      [m s-1]
      (roughness factor frc_thr_rgh_fct = 1.0 per CLM implementation)

  Step 4 — Owen's effect: saltation increases effective friction velocity:
      wnd_rfr_thr_slt = u10 * wnd_frc_thr_slt / fv      [m/s] 10-m threshold wind
      if u10 >= wnd_rfr_thr_slt:
          wnd_frc_slt = fv + 0.003 * (u10 - wnd_rfr_thr_slt)^2
      (0.003 is an empirical coefficient with implied units [s m-1])

  Step 5 — Horizontal saltation mass flux (White 1979):
      if wnd_frc_slt > wnd_frc_thr_slt:
          Q = cst_slt * forc_rho * wnd_frc_slt^3 / g
              * (1 - wnd_frc_thr_slt/wnd_frc_slt)
              * (1 + wnd_frc_thr_slt/wnd_frc_slt)^2
          Q *= fbare * mbl_bsn_fct * flx_mss_fdg_fct * liqfrac
      Units: [kg m-3][m3 s-3]/[m s-2] = [kg m-1 s-1]

  Step 6 — Vertical dust flux via sandblasting efficiency (Marticorena &
            Bergametti 1995):
      mss_frc_cly_vld = min(fclay, 0.20)
      dst_slt_flx_rat = 100 * 10^(13.4*mss_frc_cly_vld - 6)   [m-1]
                        (MB95 gives 10^(...) in cm-1; x100 converts to m-1)
      F_tot = Q * dst_slt_flx_rat                              [kg m-2 s-1]

  Step 7 — Partition into size bins via overlap of three lognormal source
            mode distributions (Bagnold, Sverdrup & Munk, Gillette 1979)
            with the four transport bins [0.1-1, 1-2.5, 2.5-5, 5-10 µm]:
      F_p[n] = sum_m( ovr_src_snk_mss[m,n] * F_tot )

Inputs (dynamic, per time step from CLM/MPAS-A):
    fv          : np.ndarray [x, y]    friction velocity [m/s]
    u10         : np.ndarray [x, y]    10-m wind speed [m/s]
    forc_rho    : np.ndarray [x, y]    near-surface air density [kg/m³]
    h2osoi_vol  : np.ndarray [x, y]    volumetric soil water content [m³/m³]
    h2osoi_liq  : np.ndarray [x, y]    liquid soil water content [kg/m²]
    h2osoi_ice  : np.ndarray [x, y]    frozen soil water content [kg/m²]
    watsat      : np.ndarray [x, y]    saturated volumetric soil water [m³/m³]
    tlai        : np.ndarray [x, y]    one-sided leaf area index [-]
    tsai        : np.ndarray [x, y]    one-sided stem area index [-]
    frac_sno    : np.ndarray [x, y]    snow cover fraction [-]

Inputs (static, loaded once at initialization):
    fclay           : np.ndarray [x, y]   clay mass fraction [-]
                                          (capped at 0.20 inside dust_emission
                                          for the sandblasting efficiency)
    mbl_bsn_fct     : np.ndarray [x, y]   basin erodibility factor [-] (= 1 in CLM4)
    is_soil         : np.ndarray [x, y]   boolean mask: soil or crop land type

Precomputed constants (from Dustini, computed once):
    tmp1            : float             threshold friction velocity pre-factor [kg^0.5 m^-0.5 s^-1]
    ovr_src_snk_mss : np.ndarray [3, 4] source-to-sink bin overlap mass fractions [-]

Output:
    F_p : np.ndarray [x, y, ndst]   surface dust emission flux [kg m^-2 s^-1]
                                    (ndst = 4 bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm)
                                    QUACS runtime converts to tendency:
                                    dx/dt = F_p / (rho_air * dz)   [kg kg-1 s-1]
                                    or volumetric source:
                                    dC/dt = F_p / dz               [kg m-3 s-1]

Note on the basin factor (mbl_bsn_fct):
    In CLM4 (DUSTMod.F90 L808-L816), mbl_bsn_fct is set to 1.0 everywhere.
    The factor originates from Ginoux et al. (2001) source function logic but
    is kept as a placeholder in DEAD/CLM4. It is retained here as an input
    to remain faithful to the Fortran and to allow future experimentation.

Note on the vegetation mask:
    CLM4 averages tlai+tsai over all PFTs at the landunit level before applying
    the VAI threshold. Here we accept a pre-averaged VAI field directly, which
    is equivalent for single-PFT use or when the user supplies the landunit mean.
"""

import numpy as np
from scipy.special import erf

# ---------------------------------------------------------------------------
# Physical constants (matching CLM4 / DUSTMod.F90)
# ---------------------------------------------------------------------------
GRAV = 9.80616  # [m s-2]  gravitational acceleration
RHO_WATER = 1000.0  # [kg m-3] density of liquid water (SHR_CONST_RHOFW)
RHO_QUARTZ = 2700.0  # [kg m-3] density of dry soil mineral
DNS_SLT = 2650.0  # [kg m-3] density of optimal saltation particles (DUSTMod L735)
DNS_AER = 2500.0  # [kg m-3] aerosol density (DUSTMod L849)

# ---------------------------------------------------------------------------
# Tuning constants (DUSTMod.F90)
# ---------------------------------------------------------------------------
CST_SLT = 2.61e0  # [-]  White (1979) saltation constant (DUSTMod L153)
FLX_MSS_FDG_FCT = 5.0e-4  # [-]  empirical mass flux tuning factor (DUSTMod L154)
VAI_MBL_THR = 0.3  # [m² m-2] VAI threshold for dust mobilization (DUSTMod L155)
FCLAY_MAX = 0.20  # [-]  clay cap for sandblasting (CTSM MassFracClay)

# ---------------------------------------------------------------------------
# Optimal saltation particle diameter (DUSTMod L734)
# ---------------------------------------------------------------------------
DMT_SLT_OPT = 75.0e-6  # [m]  optimal diameter for saltation

# ---------------------------------------------------------------------------
# Transport size bin boundaries [m]  (DUSTMod L732-L733, ndst=4)
# Bins: [0.1-1], [1-2.5], [2.5-5], [5-10]  µm
# ---------------------------------------------------------------------------
NDST = 4
DMT_GRD = np.array([0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6])  # [m] bin edges

# ---------------------------------------------------------------------------
# Source mode parameters — Kok 2011 / BSM96 (Bagnold, Sverdrup & Munk, Gillette)
# Three lognormal source modes (DUSTMod L726-L731)
# ---------------------------------------------------------------------------
DST_SRC_NBR = 3
DMT_VMA_SRC = np.array([0.832e-6, 4.82e-6, 19.38e-6])  # [m]  mass median diameters
GSD_ANL_SRC = np.array([2.10, 1.90, 1.60])  # [-]  geometric std deviations
MSS_FRC_SRC = np.array([0.036, 0.957, 0.007])  # [-]  mass fractions (sum=1)


# ---------------------------------------------------------------------------
# Precomputed constants (computed by precompute_constants(), mirrors Dustini)
# ---------------------------------------------------------------------------
def precompute_constants():
    """
    Compute the two quantities that Dustini() calculates once at model start:

      1. tmp1  — the air-density-independent prefactor of the threshold
                 friction velocity for saltation, following Marticorena &
                 Bergametti (1995) / Iversen & White (1982), evaluated at
                 the optimal saltation diameter dmt_slt_opt = 75 µm.
                 (DUSTMod L791-L806)

      2. ovr_src_snk_mss[m, n]  — mass overlap of each source lognormal
                 mode m with each transport sink bin n, used to partition
                 the total vertical dust flux into the four CESM bins.
                 (DUSTMod L777-L785)

    Returns
    -------
    tmp1            : float             [kg^0.5 m^-0.5 s^-1]
    ovr_src_snk_mss : np.ndarray [3,4]  [-]
    """
    # --- tmp1: threshold friction velocity prefactor -------------------------
    # Reynolds number approximation at optimal diameter (DUSTMod L791)
    # (100 * DMT_SLT_OPT converts m -> cm, as in the Fortran)
    ryn_nbr_frc_thr_prx_opt = 0.38 + 1331.0 * (100.0 * DMT_SLT_OPT) ** 1.56

    if ryn_nbr_frc_thr_prx_opt < 10.0:
        # Eq. for Re < 10 (DUSTMod L797-L798)
        ryn_nbr_frc_thr_opt_fnc = -1.0 + 1.928 * ryn_nbr_frc_thr_prx_opt**0.0922
        ryn_nbr_frc_thr_opt_fnc = (0.1291**2) / ryn_nbr_frc_thr_opt_fnc
    else:
        # Eq. for Re >= 10 (DUSTMod L800-L801)
        ryn_nbr_frc_thr_opt_fnc = 1.0 - 0.0858 * np.exp(-0.0617 * (ryn_nbr_frc_thr_prx_opt - 10.0))
        ryn_nbr_frc_thr_opt_fnc = (0.120**2) * (ryn_nbr_frc_thr_opt_fnc**2)

    # Inter-particle cohesive forces factor and density factor (DUSTMod L804-L806)
    icf_fct = 1.0 + 6.0e-7 / (DNS_SLT * GRAV * DMT_SLT_OPT**2.5)
    dns_fct = DNS_SLT * GRAV * DMT_SLT_OPT  # [kg m-1 s-2]
    tmp1 = np.sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)  # [kg^0.5 m^-0.5 s^-1]

    # --- ovr_src_snk_mss: lognormal source-to-sink bin overlap -------------
    # For each source mode m, compute the mass fraction landing in each
    # transport bin n using the error function integral of the lognormal.
    # (DUSTMod L777-L785)
    ovr_src_snk_mss = np.zeros((DST_SRC_NBR, NDST))
    for m in range(DST_SRC_NBR):
        sqrt2lngsd = np.sqrt(2.0) * np.log(GSD_ANL_SRC[m])
        for n in range(NDST):
            lndmax = np.log(DMT_GRD[n + 1] / DMT_VMA_SRC[m])
            lndmin = np.log(DMT_GRD[n] / DMT_VMA_SRC[m])
            ovr_src_snk_frc = 0.5 * (erf(lndmax / sqrt2lngsd) - erf(lndmin / sqrt2lngsd))
            ovr_src_snk_mss[m, n] = ovr_src_snk_frc * MSS_FRC_SRC[m]

    return tmp1, ovr_src_snk_mss


# Compute module-level constants once at import time (equivalent to Dustini)
TMP1, OVR_SRC_SNK_MSS = precompute_constants()


# ---------------------------------------------------------------------------
# Soil moisture threshold: Fécan et al. (1999) with Zender (2003a) tuning
# ---------------------------------------------------------------------------
def gwc_threshold(fclay):
    """
    Gravimetric soil water content threshold above which soil moisture
    raises the threshold friction velocity for saltation (DEAD / CLM).

    Fécan et al. (1999), with fclay as a mass fraction:
        w' = 0.17 * fclay + 0.14 * fclay^2                 [kg kg-1]
    Zender et al. (2003a) tune this by a = 1/fclay, giving the form used
    in CLM/CTSM (ThresholdSoilMoistZender2003 = 0.17 + 0.14*clay[%]*0.01):

        gwc_thr = 0.17 + 0.14 * fclay                      [kg kg-1]

    Parameters
    ----------
    fclay : np.ndarray [x, y]   soil clay mass fraction [-]

    Returns
    -------
    gwc_thr : np.ndarray [x, y]   threshold gravimetric water content [kg kg-1]
    """
    return 0.17 + 0.14 * fclay


# ---------------------------------------------------------------------------
# Core parameterization: called at every model time step
# ---------------------------------------------------------------------------
def dust_emission(
    fv,
    u10,
    forc_rho,
    h2osoi_vol,
    h2osoi_liq,
    h2osoi_ice,
    watsat,
    tlai,
    tsai,
    frac_sno,
    fclay,
    mbl_bsn_fct,
    is_soil,
    tmp1=TMP1,
    ovr_src_snk_mss=OVR_SRC_SNK_MSS,
):
    """
    Compute the DEAD (Zender et al. 2003a) dust emission flux.

    Implements the DustEmission subroutine of CTSM DustEmisZender2003.F90
    as a vectorized NumPy function over a 2-D [x, y] spatial grid.

    Physics pipeline
    ----------------
    1. Mobilizable bare-ground fraction (fbare) from VAI and snow cover.
    2. Bulk density and gravimetric soil water content.
    3. Soil-moisture threshold factor, Fécan et al. (1999) with the
       Zender (2003a) threshold gwc_thr = 0.17 + 0.14*fclay.
    4. Liquid fraction factor (inhibits emission over frozen soil).
    5. Threshold friction velocity, Marticorena & Bergametti (1995) /
       Iversen & White (1982), adjusted for moisture (and unity roughness).
    6. Owen's effect: saltation roughens the boundary layer, increasing
       the effective friction velocity above fv.
    7. Horizontal saltation mass flux, White (1979)            [kg m-1 s-1].
    8. Vertical dust flux via sandblasting efficiency, Marticorena &
       Bergametti (1995): 100*10^(13.4*min(fclay,0.2) - 6)    [m-1].
    9. Partition total flux into four transport bins via the precomputed
       lognormal overlap factors ovr_src_snk_mss.

    Parameters
    ----------
    fv           : np.ndarray [x, y]   friction velocity [m/s]
    u10          : np.ndarray [x, y]   10-m wind speed [m/s]
    forc_rho     : np.ndarray [x, y]   near-surface air density [kg/m³]
    h2osoi_vol   : np.ndarray [x, y]   volumetric soil moisture [m³/m³]
    h2osoi_liq   : np.ndarray [x, y]   liquid soil water [kg/m²]
    h2osoi_ice   : np.ndarray [x, y]   frozen soil water [kg/m²]
    watsat       : np.ndarray [x, y]   saturated volumetric soil water [m³/m³]
    tlai         : np.ndarray [x, y]   one-sided leaf area index [-]
    tsai         : np.ndarray [x, y]   one-sided stem area index [-]
    frac_sno     : np.ndarray [x, y]   snow cover fraction [-], in [0, 1]
    fclay        : np.ndarray [x, y]   clay mass fraction [-]
    mbl_bsn_fct  : np.ndarray [x, y]   basin erodibility factor [-] (= 1 in CLM4)
    is_soil      : np.ndarray [x, y]   bool, True for soil/crop land type
    tmp1         : float                threshold friction velocity prefactor
                                        (from precompute_constants)
    ovr_src_snk_mss : np.ndarray [3,4] source-to-sink mass overlap fractions
                                        (from precompute_constants)

    Returns
    -------
    F_p : np.ndarray [x, y, ndst]   dust emission flux [kg m-2 s-1]
                                     ndst = 4 bins: 0.1-1, 1-2.5, 2.5-5, 5-10 µm
    """
    nx, ny = fv.shape

    # ------------------------------------------------------------------
    # Step 1: Mobilizable bare-ground fraction (DUSTMod L311-L320)
    # fbare decreases linearly from 1→0 as VAI goes 0→vai_mbl_thr.
    # Zero over ice/lake; attenuated by snow cover.
    # ------------------------------------------------------------------
    vai = tlai + tsai  # vegetation area index
    fbare = np.where(
        is_soil,
        np.maximum(0.0, 1.0 - vai / VAI_MBL_THR) * (1.0 - frac_sno),
        0.0,
    )

    # ------------------------------------------------------------------
    # Step 2: Bulk density and gravimetric soil water content (DUSTMod L366-L367)
    # ------------------------------------------------------------------
    bd = (1.0 - watsat) * RHO_QUARTZ  # [kg m-3] bulk density of dry soil
    gwc_sfc = h2osoi_vol * RHO_WATER / bd  # [kg kg-1] gravimetric water content

    # ------------------------------------------------------------------
    # Step 3: Soil-moisture threshold factor, Fécan et al. 1999
    #         (DUSTMod L368-L372); gwc_thr = 0.17 + 0.14*fclay (CTSM)
    # ------------------------------------------------------------------
    gwc_thr = gwc_threshold(fclay)
    excess = np.maximum(gwc_sfc - gwc_thr, 0.0)
    frc_thr_wet_fct = np.where(
        gwc_sfc > gwc_thr,
        np.sqrt(1.0 + 1.21 * (100.0 * excess) ** 0.68),
        1.0,
    )

    # ------------------------------------------------------------------
    # Step 4: Liquid fraction — inhibits emission over frozen soil
    #         (DUSTMod L376)
    # ------------------------------------------------------------------
    liqfrac = np.clip(
        h2osoi_liq / (h2osoi_ice + h2osoi_liq + 1.0e-6),
        0.0,
        1.0,
    )

    # ------------------------------------------------------------------
    # Step 5: Threshold friction velocity for saltation
    #         Iversen & White (1982) / Marticorena & Bergametti (1995)
    #         evaluated at optimal saltation diameter via the pre-factor
    #         tmp1, then adjusted for moisture (DUSTMod L384)
    #         frc_thr_rgh_fct = 1.0 (no roughness correction in CLM4)
    #         [kg^0.5 m^-0.5 s^-1] / [kg^0.5 m^-1.5] = [m s-1]
    # ------------------------------------------------------------------
    wnd_frc_thr_slt = tmp1 / np.sqrt(forc_rho) * frc_thr_wet_fct  # [m/s]

    # ------------------------------------------------------------------
    # Step 6: Owen's effect — saltation increases effective u* (DUSTMod L395-L405)
    #         10-m reference wind threshold:
    #             wnd_rfr_thr_slt = u10 * wnd_frc_thr_slt / fv
    #         Saltation-enhanced friction velocity:
    #             wnd_frc_slt = fv + 0.003 * (u10 - wnd_rfr_thr_slt)^2
    # ------------------------------------------------------------------
    # Avoid division by zero for calm cells
    fv_safe = np.maximum(fv, 1.0e-6)
    wnd_rfr_thr_slt = u10 * wnd_frc_thr_slt / fv_safe  # [m/s]
    wnd_rfr_dlt = np.maximum(u10 - wnd_rfr_thr_slt, 0.0)  # [m/s]
    wnd_frc_slt_dlt = 0.003 * wnd_rfr_dlt**2  # [m/s]
    wnd_frc_slt = np.where(
        u10 >= wnd_rfr_thr_slt,
        fv + wnd_frc_slt_dlt,
        fv,
    )

    # ------------------------------------------------------------------
    # Step 7: Horizontal saltation mass flux, White (1979) (DUSTMod L410-L421)
    #
    #   Q = C_slt * rho * u*_slt^3 / g
    #       * (1 - u*_thr/u*_slt) * (1 + u*_thr/u*_slt)^2
    #   Units: [kg m-3][m3 s-3]/[m s-2] = [kg m-1 s-1]
    #
    #   Applied only where wnd_frc_slt > wnd_frc_thr_slt,
    #   then multiplied by fbare, basin factor, tuning factor, liqfrac.
    # ------------------------------------------------------------------
    emit_mask = fbare > 0.0

    wnd_frc_rat = np.where(
        emit_mask & (wnd_frc_slt > wnd_frc_thr_slt),
        wnd_frc_thr_slt / wnd_frc_slt,
        1.0,  # ratio = 1 → Q = 0
    )

    Q = np.where(
        emit_mask & (wnd_frc_slt > wnd_frc_thr_slt),
        (
            CST_SLT
            * forc_rho
            * wnd_frc_slt**3
            / GRAV
            * (1.0 - wnd_frc_rat)
            * (1.0 + wnd_frc_rat) ** 2
            * fbare
            * mbl_bsn_fct
            * FLX_MSS_FDG_FCT
            * liqfrac
        ),
        0.0,
    )  # [kg m-1 s-1]

    # ------------------------------------------------------------------
    # Step 8: Vertical dust flux via sandblasting efficiency
    #         Marticorena & Bergametti (1995) (DUSTMod L428-L429)
    #
    #   mss_frc_cly_vld = min(fclay, 0.20)                (CTSM MassFracClay)
    #   dst_slt_flx_rat = 100 * 10^(13.4*mss_frc_cly_vld - 6)    [m-1]
    #   F_tot = Q * dst_slt_flx_rat                               [kg m-2 s-1]
    # ------------------------------------------------------------------
    mss_frc_cly_vld = np.minimum(fclay, FCLAY_MAX)
    dst_slt_flx_rat = 100.0 * 10.0 ** (13.4 * mss_frc_cly_vld - 6.0)  # [m-1]
    F_tot = Q * dst_slt_flx_rat  # [kg m-2 s-1]

    # ------------------------------------------------------------------
    # Step 9: Partition total flux into four transport bins (DUSTMod L438-L447)
    #         F_p[n] = sum_m( ovr_src_snk_mss[m, n] * F_tot )
    # ------------------------------------------------------------------
    F_p = np.zeros((nx, ny, NDST))
    for n in range(NDST):
        bin_frac = ovr_src_snk_mss[:, n].sum()  # sum over source modes m
        F_p[:, :, n] = F_tot * bin_frac

    return F_p  # [kg m-2 s-1], shape [x, y, ndst=4]
