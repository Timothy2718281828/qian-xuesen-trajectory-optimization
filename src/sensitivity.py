"""
sensitivity.py -- M7: Sensitivity analysis of optimal trajectory parameters

Analyses how Dv_total responds to perturbations in the key design variables:
  - rp (perihelion distance)
  - rm (Moon flyby distance)
  - t0 (launch day)
  - theta_0 (Moon initial phase uncertainty)

Produces numerical partial derivatives and generates data for plotting.

Usage:
    python src/sensitivity.py
"""

import numpy as np
import json
import os
import sys

from physical_constants import (
    AU, R_MOON, R_SUN,
    C_RM_MIN, C_RM_MAX, C_RP_MIN, C_RP_MAX,
    SECONDS_PER_DAY,
)
from single_trajectory import (
    solve_best_trajectory,
    solve_trajectory,
    moon_phase_at_day,
    THETA_MOON_J2026,
    N_MOON_PER_DAY,
)
from lunar_flyby import direct_launch_delta_v, design_trajectory

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


def partial_dv_rp(t0_day, rm, rp_center, delta_rp=None, n_points=21):
    """
    d(Dv_total)/d(rp) curve around a centre point.

    Returns arrays (rp_values_AU, dv_values, dv_direct_values).
    """
    if delta_rp is None:
        delta_rp = 0.15 * AU

    rp_lo = max(rp_center - delta_rp, C_RP_MIN)
    rp_hi = min(rp_center + delta_rp, C_RP_MAX)
    rp_vals = np.linspace(rp_lo, rp_hi, n_points)

    rp_au, dv_a, dv_d = [], [], []
    for rp in rp_vals:
        b = solve_best_trajectory(t0_day, rm, rp)
        d = direct_launch_delta_v(rp)
        if b is not None:
            rp_au.append(rp / AU)
            dv_a.append(b['delta_v_total'])
            dv_d.append(d['delta_v_total'])

    return np.array(rp_au), np.array(dv_a), np.array(dv_d)


def partial_dv_rm(t0_day, rm_center, rp, n_points=21):
    """
    d(Dv_total)/d(rm) curve.

    Returns arrays (rm_values_km, dv_values).
    """
    rm_lo = max(C_RM_MIN, rm_center * 0.3)
    rm_hi = min(C_RM_MAX, rm_center * 5.0)
    rm_vals = np.logspace(np.log10(rm_lo), np.log10(rm_hi), n_points)

    rm_out, dv_out = [], []
    for rm in rm_vals:
        b = solve_best_trajectory(t0_day, rm, rp)
        if b is not None:
            rm_out.append(rm)
            dv_out.append(b['delta_v_total'])

    return np.array(rm_out), np.array(dv_out)


def partial_dv_t0(rm, rp, t0_center, window=14, step=1):
    """
    d(Dv_total)/d(t0) curve around a launch day.

    Returns arrays (day_offsets, dv_values).
    """
    days = np.arange(t0_center - window, t0_center + window + 1, step)
    days = days[(days >= 0) & (days < 365)]

    day_out, dv_out = [], []
    for day in days:
        b = solve_best_trajectory(int(day), rm, rp)
        if b is not None:
            day_out.append(int(day))
            dv_out.append(b['delta_v_total'])

    return np.array(day_out), np.array(dv_out)


def partial_dv_theta(t0_day, rm, rp, delta_deg=10.0, n_points=21):
    """
    d(Dv_total)/d(theta0) -- sensitivity to Moon initial phase uncertainty.

    Directly varies the Moon phase fed to design_trajectory.
    """
    base_phase = moon_phase_at_day(t0_day)
    d_theta = np.linspace(-np.radians(delta_deg), np.radians(delta_deg), n_points)

    theta_deg, dv_out = [], []
    for dt in d_theta:
        phase = (base_phase + dt) % (2.0 * np.pi)
        for side in ('leading', 'trailing'):
            r = design_trajectory(rp, rm, side, phase)
            if r is not None:
                theta_deg.append(np.degrees(dt))
                dv_out.append(r['delta_v_total'])
                break

    return np.array(theta_deg), np.array(dv_out)


def numerical_gradient(x, y):
    """Central-difference gradient dy/dx, with forward/backward at endpoints."""
    if len(x) < 2:
        return np.array([0.0])
    grad = np.zeros_like(y)
    for i in range(len(x)):
        if i == 0:
            grad[i] = (y[1] - y[0]) / (x[1] - x[0])
        elif i == len(x) - 1:
            grad[i] = (y[-1] - y[-2]) / (x[-1] - x[-2])
        else:
            grad[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    return grad


def full_sensitivity_report(t0_day, rm, rp, save=True):
    """
    Compute all partial-derivative curves and summary statistics.

    Returns a dict with curves and key metrics.
    """
    print("=" * 60)
    print(f"  M7: Sensitivity Analysis  (day={t0_day}, rm={rm:.0f}, rp={rp/AU:.3f} AU)")
    print("=" * 60)

    report = {'params': {'t0_day': t0_day, 'rm': rm, 'rp_AU': rp / AU}}

    # 1. rp sensitivity
    print("\n  [1/4] Dv vs rp ...")
    rp_au, dv_rp, dv_dir = partial_dv_rp(t0_day, rm, rp, n_points=25)
    if len(rp_au) >= 2:
        grad_rp = numerical_gradient(rp_au, dv_rp)
        report['rp'] = {
            'rp_AU': rp_au.tolist(), 'dv': dv_rp.tolist(),
            'dv_direct': dv_dir.tolist(), 'grad_dv_drp': grad_rp.tolist(),
        }
        idx_mid = len(rp_au) // 2
        print(f"    d(Dv)/d(rp) ≈ {grad_rp[idx_mid]:.2f} km/s per AU  (at rp={rp_au[idx_mid]:.3f} AU)")
    else:
        report['rp'] = None
        print("    Not enough data points.")

    # 2. rm sensitivity
    print("\n  [2/4] Dv vs rm ...")
    rm_arr, dv_rm = partial_dv_rm(t0_day, rm, rp, n_points=25)
    if len(rm_arr) >= 2:
        grad_rm = numerical_gradient(rm_arr, dv_rm)
        report['rm'] = {
            'rm_km': rm_arr.tolist(), 'dv': dv_rm.tolist(),
            'grad_dv_drm': grad_rm.tolist(),
        }
        idx_mid = len(rm_arr) // 2
        print(f"    d(Dv)/d(rm) ≈ {grad_rm[idx_mid]*1000:.4f} km/s per 1000 km  (at rm={rm_arr[idx_mid]:.0f} km)")
    else:
        report['rm'] = None
        print("    Not enough data points.")

    # 3. t0 sensitivity
    print("\n  [3/4] Dv vs t0 ...")
    days, dv_t0 = partial_dv_t0(rm, rp, t0_day, window=14)
    if len(days) >= 2:
        grad_t0 = numerical_gradient(days.astype(float), dv_t0)
        report['t0'] = {
            'days': days.tolist(), 'dv': dv_t0.tolist(),
            'grad_dv_dt0': grad_t0.tolist(),
        }
        idx_mid = len(days) // 2
        print(f"    d(Dv)/d(t0) ≈ {grad_t0[idx_mid]:.4f} km/s per day  (at day={days[idx_mid]})")
    else:
        report['t0'] = None
        print("    Not enough data points.")

    # 4. Moon phase sensitivity
    print("\n  [4/4] Dv vs theta_0 ...")
    theta_deg, dv_theta = partial_dv_theta(t0_day, rm, rp, delta_deg=15.0, n_points=31)
    if len(theta_deg) >= 2:
        grad_theta = numerical_gradient(theta_deg, dv_theta)
        report['theta'] = {
            'delta_deg': theta_deg.tolist(), 'dv': dv_theta.tolist(),
            'grad_dv_dtheta': grad_theta.tolist(),
        }
        idx_mid = len(theta_deg) // 2
        print(f"    d(Dv)/d(theta) ≈ {grad_theta[idx_mid]:.4f} km/s per deg  (at delta_theta=0)")
    else:
        report['theta'] = None
        print("    Not enough data points.")

    # Summary
    print("\n  === Summary ===")
    dv_base = solve_best_trajectory(t0_day, rm, rp)
    direct = direct_launch_delta_v(rp)
    if dv_base is not None:
        sv = direct['delta_v_total'] - dv_base['delta_v_total']
        print(f"  Baseline Dv_total  = {dv_base['delta_v_total']:.4f} km/s")
        print(f"  Direct   Dv_total  = {direct['delta_v_total']:.4f} km/s")
        print(f"  Savings            = {sv:.4f} km/s ({sv/direct['delta_v_total']*100:.2f}%)")

    # Assess practical significance of Moon phase uncertainty
    if report.get('theta') and len(report['theta']['dv']) > 1:
        dv_arr = np.array(report['theta']['dv'])
        dv_range = dv_arr.max() - dv_arr.min()
        print(f"  Dv range over ±15° phase: {dv_range:.4f} km/s")
        if sv > 0:
            print(f"  Phase error > ±{abs(sv / grad_theta[idx_mid]):.1f}° would negate savings" if abs(grad_theta[idx_mid]) > 1e-8 else "  Phase sensitivity negligible")

    report['summary'] = {
        'dv_base': dv_base['delta_v_total'] if dv_base else None,
        'dv_direct': direct['delta_v_total'],
        'savings': sv if dv_base else 0.0,
    }

    if save:
        path = os.path.join(DATA_DIR, 'sensitivity.json')
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        print(f"\n  Saved: {path}")

    return report


if __name__ == '__main__':
    # Use optimal parameters from C3 scan
    optimal_path = os.path.join(DATA_DIR, 'scan_2026_optimal.json')
    if os.path.exists(optimal_path):
        with open(optimal_path, 'r') as f:
            opt = json.load(f)
        t0 = opt['t0_day']
        rm = opt['rm_optimal']
        rp = opt['rp_optimal']
    else:
        t0, rm, rp = 52, 5000.0, 0.2 * AU

    report = full_sensitivity_report(t0, rm, rp)

    # Also run for rp=0.2 AU (classic Qian Xuesen case)
    if rp / AU > 0.25:
        print("\n\n" + "=" * 60)
        print("  Sensitivity for classic case (rp=0.2 AU)")
        print("=" * 60)
        full_sensitivity_report(t0, rm, 0.2 * AU, save=False)
