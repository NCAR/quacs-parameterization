"""Step 1 driver: diagnose the pyroCb flag."""

from __future__ import annotations

import math

import numpy as np

from quacs.plumerise.tables import lookup_heat_flux

MIN_FIRE_SIZE = 0.0
F_CONV = 0.55
MOISTURE_RATIO = 15.0
PFT_CONST = 397.3
TEMPERATURE_BUFFER = 0.5
SP_BETA_MAX = 0.25
DELTA_BETA = 0.01
MAX_PLUME_TOP_T_C = -20.0

R_DRY_AIR = 287.04
R_WATER_VAPOR = 461.55
CP_DRY_AIR = 1004.0
EPSILON = R_DRY_AIR / R_WATER_VAPOR


def _validate_profile(z, p, t, u, v, qv):
    sizes = {
        len(z),
        len(p),
        len(t),
        len(u),
        len(v),
        len(qv),
    }
    if len(sizes) != 1 or next(iter(sizes)) < 2:
        raise ValueError("meteorological profile arrays must share length >= 2")
    if np.any(np.diff(z) <= 0.0):
        raise ValueError("z must be strictly increasing")
    if np.any(np.diff(p) >= 0.0):
        raise ValueError("p must decrease with height")
    if np.any(p <= 0.0):
        raise ValueError("p must be positive")


def _c_to_k(t_c):
    return t_c + 273.15


def _k_to_c(t_k):
    return t_k - 273.15


def _theta_kelvin(pressure_hpa, t_c):
    return _c_to_k(t_c) * (1000.0 / pressure_hpa) ** (R_DRY_AIR / CP_DRY_AIR)


def _temperature_c_from_theta(theta_kelvin, pressure_hpa):
    return _k_to_c(theta_kelvin * (pressure_hpa / 1000.0) ** (R_DRY_AIR / CP_DRY_AIR))


def _vapor_pressure_liquid_water(t_c):
    return 6.1037 * math.exp(17.641 * t_c / (t_c + 243.27))


def _dew_point_from_vapor_pressure(vp_hpa):
    if vp_hpa <= 0.0:
        raise ValueError("vapor pressure must be positive")
    a = math.log(vp_hpa / 6.1037) / 17.641
    return a * 243.27 / (1.0 - a)


def _specific_humidity_from_dewpoint(dp_c, pressure_hpa):
    return _vapor_pressure_liquid_water(dp_c) / pressure_hpa * EPSILON


def _dew_point_from_p_and_specific_humidity(pressure_hpa, q):
    vp_hpa = max(q, 1.0e-12) * pressure_hpa / EPSILON
    return _dew_point_from_vapor_pressure(vp_hpa)


def _latent_heat_of_condensation(t_c):
    if t_c < -100.0 or t_c > 60.0:
        raise ValueError("temperature outside latent-heat approximation range")
    return (2500.8 - 2.36 * t_c + 0.0016 * t_c**2 - 0.00006 * t_c**3) * 1000.0


def _temperature_kelvin_at_lcl(t_c, dp_c):
    if dp_c >= t_c:
        return _c_to_k(t_c)
    celsius_lcl = dp_c - (0.001296 * dp_c + 0.1963) * (t_c - dp_c)
    return _c_to_k(celsius_lcl)


def _press_and_temp_k_at_lcl(t_c, dp_c, pressure_hpa):
    if dp_c >= t_c:
        return pressure_hpa, _c_to_k(t_c)
    t_lcl = _temperature_kelvin_at_lcl(t_c, dp_c)
    p_lcl = pressure_hpa * (t_lcl / _c_to_k(t_c)) ** (CP_DRY_AIR / R_DRY_AIR)
    return p_lcl, t_lcl


def _theta_e_kelvin(t_c, dp_c, pressure_hpa):
    latent_heat = _latent_heat_of_condensation(t_c)
    theta = _theta_kelvin(pressure_hpa, t_c)
    t_lcl = _temperature_kelvin_at_lcl(t_c, dp_c)
    q = _specific_humidity_from_dewpoint(dp_c, pressure_hpa)
    return theta * (1.0 + latent_heat * q / (CP_DRY_AIR * t_lcl))


def _theta_e_saturated_kelvin(pressure_hpa, t_c):
    latent_heat = _latent_heat_of_condensation(t_c)
    theta = _theta_kelvin(pressure_hpa, t_c)
    q_sat = _specific_humidity_from_dewpoint(t_c, pressure_hpa)
    return theta * (1.0 + latent_heat * q_sat / (CP_DRY_AIR * _c_to_k(t_c)))


def _find_root(func, low_val, high_val):
    if low_val > high_val:
        low_val, high_val = high_val, low_val

    f_low = func(low_val)
    f_high = func(high_val)
    if f_low is None or f_high is None:
        return None
    if f_high * f_low > 0.0:
        return None

    mid_val = (high_val - low_val) / 2.0 + low_val
    for _ in range(50):
        f_mid = func(mid_val)
        if f_mid is None:
            return None
        if f_mid * f_low > 0.0:
            low_val = mid_val
            f_low = f_mid
        else:
            high_val = mid_val
        if abs(high_val - low_val) < 1.0e-10:
            break
        mid_val = (high_val - low_val) / 2.0 + low_val
    return mid_val


def _temperature_c_from_theta_e_saturated_and_pressure(pressure_hpa, theta_e):
    def diff(t_c):
        try:
            return _theta_e_saturated_kelvin(pressure_hpa, t_c) - theta_e
        except ValueError:
            return None

    return _find_root(diff, -80.0, 50.0)


def _virtual_temperature_c(t_c, dp_c, pressure_hpa):
    mixing_ratio = EPSILON * _vapor_pressure_liquid_water(dp_c) / (
        pressure_hpa - _vapor_pressure_liquid_water(dp_c)
    )
    theta = _theta_kelvin(pressure_hpa, t_c)
    virtual_theta = theta * (1.0 + mixing_ratio / EPSILON) / (1.0 + mixing_ratio)
    return _temperature_c_from_theta(virtual_theta, pressure_hpa)


def _interp(z, values, z_target):
    return float(np.interp(z_target, z, values))


def _theta_q_to_p_t(theta, q):
    def diff(pressure_hpa):
        t_c = _temperature_c_from_theta(theta, pressure_hpa)
        dp_c = _dew_point_from_p_and_specific_humidity(pressure_hpa, q)
        return t_c - dp_c

    pressure_hpa = _find_root(diff, 1080.0, 100.0)
    if pressure_hpa is None:
        return None, None
    return pressure_hpa, _temperature_c_from_theta(theta, pressure_hpa)


def _entrained_mixed_layer(pressure, z, t_c, qv):
    p_sfc = float(pressure[0])
    max_pressure = p_sfc - 50.0
    z_rel = z - z[0]

    rows = [
        (p, height, _theta_kelvin(p, temp), q)
        for p, height, temp, q in zip(pressure, z_rel, t_c, qv)
        if np.isfinite(p) and np.isfinite(height) and np.isfinite(temp) and np.isfinite(q)
    ]
    if len(rows) < 3:
        raise ValueError("not enough valid levels for entrained mixed layer")

    sum_theta = 0.0
    sum_q = 0.0
    avg_rows = []
    for row0, row1 in zip(rows, rows[1:]):
        _p0, h0, theta0, q0 = row0
        p1, h1, theta1, q1 = row1
        if h1 <= h0:
            raise ValueError("z must increase through the mixed layer")

        dh = h1 - h0
        sum_theta += (theta0 * h0 + theta1 * h1) * dh
        sum_q += (q0 * h0 + q1 * h1) * dh
        if p1 >= max_pressure:
            continue
        h_sq = h1 * h1
        avg_theta = sum_theta / h_sq
        avg_q = sum_q / h_sq
        t_mixed_c = _temperature_c_from_theta(avg_theta, p_sfc)
        dp_mixed_c = _dew_point_from_p_and_specific_humidity(p_sfc, avg_q)
        p_lcl, _t_lcl = _press_and_temp_k_at_lcl(t_mixed_c, dp_mixed_c, p_sfc)
        avg_rows.append((p1, h1, avg_theta, avg_q, p_lcl))

    if len(avg_rows) < 2:
        raise ValueError("could not define a 50 hPa entrained mixed layer")

    for row0, row1 in zip(avg_rows, avg_rows[1:]):
        p0, _h0, avg_theta0, avg_q0, p_lcl0 = row0
        p1, _h1, _avg_theta1, _avg_q1, p_lcl1 = row1
        if p0 > p_lcl0 and p1 < p_lcl1:
            return avg_theta0, avg_q0, p_sfc

    raise ValueError("could not find entrained mixed-layer LCL crossing")


def _min_temperature_diff_to_max_cloud_top_temperature(pressure, t_c, qv, starting_pressure, theta_e):
    diffs = []
    selected_count = 0
    for p, temp, q in zip(pressure, t_c, qv):
        if p > starting_pressure:
            continue
        parcel_t = _temperature_c_from_theta_e_saturated_and_pressure(p, theta_e)
        if parcel_t is None:
            continue
        if parcel_t < MAX_PLUME_TOP_T_C and selected_count >= 2:
            break
        dp_env = _dew_point_from_p_and_specific_humidity(p, q)
        env_virtual_t = _virtual_temperature_c(temp, dp_env, p)
        parcel_virtual_t = _virtual_temperature_c(parcel_t, parcel_t, p)
        diffs.append(parcel_virtual_t - env_virtual_t)
        selected_count += 1

    if not diffs:
        return None
    return min(diffs)


def _is_free_convecting(pressure, t_c, qv, starting_pressure, theta_e):
    min_diff = _min_temperature_diff_to_max_cloud_top_temperature(
        pressure,
        t_c,
        qv,
        starting_pressure,
        theta_e,
    )
    if min_diff is None:
        return False
    return min_diff >= TEMPERATURE_BUFFER


def _free_convection_level(pressure, z, t_c, qv, theta_ml, q_ml):
    def apply_beta(beta):
        theta_sp = (1.0 + beta) * theta_ml
        q_sp = q_ml + beta / MOISTURE_RATIO / 1000.0 * theta_ml
        return theta_sp, q_sp

    low_beta = math.nan
    high_beta = math.nan
    beta = 0.0
    while beta <= SP_BETA_MAX:
        theta_sp, q_sp = apply_beta(beta)
        p_sp, t_sp = _theta_q_to_p_t(theta_sp, q_sp)
        if p_sp is None or t_sp is None:
            beta += DELTA_BETA
            continue

        theta_e = _theta_e_kelvin(t_sp, t_sp, p_sp)
        free_convecting = _is_free_convecting(pressure, t_c, qv, p_sp, theta_e)
        if not free_convecting and math.isnan(high_beta):
            low_beta = beta
        elif free_convecting and math.isnan(high_beta):
            high_beta = beta
            break
        beta += DELTA_BETA

    if math.isnan(low_beta) or math.isnan(high_beta):
        raise ValueError("could not bracket free-convection beta")

    def root_func(beta_value):
        theta_sp, q_sp = apply_beta(beta_value)
        p_sp, t_sp = _theta_q_to_p_t(theta_sp, q_sp)
        if p_sp is None or t_sp is None:
            return None
        theta_e = _theta_e_kelvin(t_sp, t_sp, p_sp)
        min_diff = _min_temperature_diff_to_max_cloud_top_temperature(
            pressure,
            t_c,
            qv,
            p_sp,
            theta_e,
        )
        if min_diff is None:
            return None
        return min_diff - TEMPERATURE_BUFFER

    beta_fc = _find_root(root_func, low_beta, high_beta)
    if beta_fc is None:
        raise ValueError("could not solve free-convection beta")

    theta_fc, q_sp = apply_beta(beta_fc)
    dtheta_fc = theta_fc - theta_ml
    p_fc, _t_fc = _theta_q_to_p_t(theta_fc, q_sp)
    if p_fc is None:
        raise ValueError("could not solve free-convection pressure")

    high_idx = None
    low_idx = None
    for idx, p in enumerate(pressure):
        if p > p_fc:
            low_idx = idx
        elif p <= p_fc:
            high_idx = idx
            break
    if low_idx is None or high_idx is None:
        raise ValueError("free-convection pressure is outside profile")

    z_fc = _interp(
        pressure[high_idx : low_idx - 1 : -1],
        z[high_idx : low_idx - 1 : -1],
        p_fc,
    )
    return z_fc - z[0], p_fc, theta_fc, dtheta_fc


def _entrained_layer_mean_wind_speed(z_fc, pressure, z, u, v):
    z_top = z[0] + z_fc
    rows = [
        (p, wind_u, wind_v)
        for p, height, wind_u, wind_v in zip(pressure, z, u, v)
        if height <= z_top and np.isfinite(p) and np.isfinite(wind_u) and np.isfinite(wind_v)
    ]
    if len(rows) < 2:
        raise ValueError("not enough wind levels below free-convection height")

    sum_p = 0.0
    sum_u = 0.0
    sum_v = 0.0
    for row0, row1 in zip(rows, rows[1:]):
        p0, u0, v0 = row0
        p1, u1, v1 = row1
        dp = p1 - p0
        sum_u += (u0 * p0 + u1 * p1) * dp
        sum_v += (v0 * p0 + v1 * p1) * dp
        sum_p += (p0 + p1) * dp

    if sum_p == 0.0:
        raise ValueError("cannot compute pressure-weighted wind speed")
    avg_u = sum_u / sum_p
    avg_v = sum_v / sum_p
    return float(math.sqrt(avg_u**2 + avg_v**2))


def _pft_formula(z_fc, p_fc, mean_wind, dtheta_fc, theta_fc, p_sfc):
    z_fc_km = z_fc / 1000.0
    p_c = p_sfc - (p_sfc - p_fc) / (1.0 + 0.32 * 0.4)
    density = p_c / 10.0 / (R_DRY_AIR * theta_fc) * (
        (1000.0 / p_c) ** (R_DRY_AIR / CP_DRY_AIR)
    )
    return float(PFT_CONST * density * (z_fc_km**2) * mean_wind * dtheta_fc)


def compute_PFT(z, p, t, u, v, qv):
    _validate_profile(z, p, t, u, v, qv)
    pressure = p
    t_c = t - 273.15
    theta_ml, q_ml, p_sfc = _entrained_mixed_layer(pressure, z, t_c, qv)
    z_fc, p_fc, theta_fc, dtheta_fc = _free_convection_level(
        pressure,
        z,
        t_c,
        qv,
        theta_ml,
        q_ml,
    )
    mean_wind = _entrained_layer_mean_wind_speed(z_fc, pressure, z, u, v)
    return _pft_formula(z_fc, p_fc, mean_wind, dtheta_fc, theta_fc, p_sfc)


def compute_qplume(fire_size_mean, fire_size_std, vegetation_class):
    fire_size_high = max(fire_size_mean + fire_size_std, 0.0)
    heat_flux, _source_key = lookup_heat_flux(vegetation_class)
    heat_kw = fire_size_high * heat_flux.heat_flux_kw_m2
    return float(1.0e-6 * F_CONV * heat_kw)


def diagnose_pyrocb_flag(PFT, qplume):
    return qplume >= PFT


def pyrocb_flag_driver(
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
    if fire_size_mean <= MIN_FIRE_SIZE:
        return False

    PFT = compute_PFT(z, p, t, u, v, qv)
    qplume = compute_qplume(
        fire_size_mean,
        fire_size_std,
        vegetation_class,
    )
    pyrocb_flag = diagnose_pyrocb_flag(PFT, qplume)
    return pyrocb_flag
