"""
optimize_launch.py -- M6: Full-year optimal launch window scanner

Three-phase optimization:
  1. Coarse scan: 365 days x coarse rp grid (fixed rm)
  2. Fine scan:  Top-N days x fine (rm, rp) grid
  3. Refinement: Best candidate localized search

Uses single_trajectory.py interfaces from Claude Opus (C1/C2).

Outputs:
  - data/scan_2026_coarse.json  -- Phase-1 results
  - data/scan_2026_fine.json    -- Phase-2 results
  - data/scan_2026_optimal.json -- Final best solution

Usage:
  python src/optimize_launch.py          # Full scan (coarse only, fast)
  python src/optimize_launch.py --fine   # Coarse + fine + refine
  python src/optimize_launch.py --quick  # Quick 28-day demo scan
"""

import numpy as np
import json
import os
import sys
import time

from physical_constants import (
    AU, R_MOON, R_SUN, R_EARTH_ORBIT,
    C_RM_MIN, C_RM_MAX, C_RP_MIN, C_RP_MAX,
    SECONDS_PER_DAY,
)

# These are imported from Opus modules
from single_trajectory import (
    solve_trajectory,
    solve_best_trajectory,
    solve_or_direct,
    scan_launch_days,
    scan_rm_rp_grid,
    check_all_constraints,
    run_nbody_verification,
)
from lunar_flyby import direct_launch_delta_v

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


# ============================================================
#  Phase 1: Coarse 365-day scan
# ============================================================

def coarse_scan(n_days=365, rm_default=5000.0,
                rp_values=None, save=True):
    """
    Scan all launch days with default rm and a coarse rp grid.

    Parameters
    ----------
    n_days : int
        Number of days to scan (365 = full year)
    rm_default : float
        Default Moon flyby distance [km]
    rp_values : list or None
        Perihelion distances to try [km]. None = auto-generate.
    save : bool
        Save results to JSON

    Returns
    -------
    list of dict: [{'day', 'rm', 'rp', 'side', 'dv_total', 'savings', ...}, ...]
    """
    if rp_values is None:
        # Auto-generate: log-spaced from 0.05 AU to 0.4 AU
        rp_au = np.logspace(np.log10(0.05), np.log10(0.4), 8)
        rp_values = [r * AU for r in rp_au]

    all_results = []
    best_per_day = {}  # day -> best result

    print(f"[Coarse Scan] {n_days} days x {len(rp_values)} rp values")
    print(f"  rm fixed at {rm_default:.0f} km")
    t_start = time.time()

    for day in range(n_days):
        for rp in rp_values:
            if rp < C_RP_MIN or rp > C_RP_MAX:
                continue

            # Try both sides
            best = solve_best_trajectory(day, rm_default, rp)
            if best is not None:
                direct = direct_launch_delta_v(rp)
                savings = direct['delta_v_total'] - best['delta_v_total']
                record = {
                    'day': day,
                    'rm': rm_default,
                    'rp_target': rp,
                    'rp_target_AU': rp / AU,
                    'side': best.get('side', 'N/A'),
                    'dv_launch': best['delta_v_launch'],
                    'dv_reentry': best['delta_v_reentry'],
                    'dv_total': best['delta_v_total'],
                    'savings': savings,
                    'savings_pct': 100.0 * savings / direct['delta_v_total'],
                    'dv_direct': direct['delta_v_total'],
                    't_total_days': best.get('t_total', 0) / SECONDS_PER_DAY,
                    'rp_achieved_AU': best.get('rp_achieved', rp) / AU,
                }
                all_results.append(record)

                # Track best per day
                if day not in best_per_day or savings > best_per_day[day]['savings']:
                    best_per_day[day] = record

        if (day + 1) % 30 == 0:
            t_elapsed = time.time() - t_start
            n_results = len(all_results)
            print(f"  Day {day+1}/{n_days}: {n_results} feasible so far "
                  f"({t_elapsed:.0f}s)")

    t_total = time.time() - t_start
    print(f"\n[Coarse Scan] Complete: {len(all_results)} feasible points "
          f"in {t_total:.0f}s")

    # Find global best
    if all_results:
        best_global = max(all_results, key=lambda x: x['savings'])
        print(f"  Best: day={best_global['day']}, "
              f"rp={best_global['rp_target_AU']:.3f} AU, "
              f"dv={best_global['dv_total']:.3f} km/s, "
              f"savings={best_global['savings']:.3f} km/s "
              f"({best_global['savings_pct']:.2f}%)")

    output = {
        'config': {
            'n_days': n_days,
            'rm_default': rm_default,
            'rp_values_AU': [r / AU for r in rp_values],
        },
        'results': all_results,
        'best_per_day': {str(k): v for k, v in best_per_day.items()},
        'global_best': max(all_results, key=lambda x: x['savings']) if all_results else None,
        'n_feasible': len(all_results),
        'n_days_with_solution': len(best_per_day),
    }

    if save:
        path = os.path.join(DATA_DIR, 'scan_2026_coarse.json')
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(output, f, indent=2, default=str)
        print(f"  Saved: {path}")

    return output


# ============================================================
#  Phase 2: Fine grid scan around top candidates
# ============================================================

def fine_scan(coarse_results, top_n=10,
               rm_list=None, rp_fine_count=10, save=True):
    """
    Fine (rm, rp) grid around the best days from coarse scan.

    Parameters
    ----------
    coarse_results : dict
        Output from coarse_scan()
    top_n : int
        Number of top days to fine-scan
    rm_list : list or None
        rm values to try [km]. None = auto-generate.
    rp_fine_count : int
        Number of rp values per day

    Returns
    -------
    list of dict: Full trajectory solutions
    """
    if 'results' not in coarse_results or not coarse_results['results']:
        print("[Fine Scan] No coarse results to refine.")
        return []

    # Get top N days by savings
    all_r = coarse_results['results']
    days_seen = {}
    for r in sorted(all_r, key=lambda x: x['savings'], reverse=True):
        d = r['day']
        if d not in days_seen:
            days_seen[d] = r

    top_days = list(days_seen.keys())[:top_n]
    print(f"[Fine Scan] Top {len(top_days)} days: {top_days}")

    if rm_list is None:
        # 10 values from 2000 to 50000
        rm_list = np.logspace(np.log10(2000), np.log10(C_RM_MAX), 10)

    # For each top day, find the best rp from coarse and scan around it
    fine_results = []

    for day in top_days:
        coarse_best = days_seen[day]
        rp_center = coarse_best['rp_target']
        rp_half_range = rp_center * 0.3  # +/- 30% around best rp

        rp_fine = np.linspace(
            max(rp_center - rp_half_range, C_RP_MIN),
            min(rp_center + rp_half_range, C_RP_MAX),
            rp_fine_count
        )

        print(f"  Day {day}: rm={len(rm_list)} values, rp={rp_fine_count} values")

        for rm in rm_list:
            for rp in rp_fine:
                best = solve_best_trajectory(day, rm, rp)
                if best is not None:
                    direct = direct_launch_delta_v(rp)
                    savings = direct['delta_v_total'] - best['delta_v_total']
                    fine_results.append({
                        'day': day,
                        'rm': rm,
                        'rp_target': rp,
                        'rp_target_AU': rp / AU,
                        'side': best.get('side', 'N/A'),
                        'dv_total': best['delta_v_total'],
                        'savings': savings,
                        'savings_pct': 100.0 * savings / direct['delta_v_total'],
                        'dv_direct': direct['delta_v_total'],
                        'result': best,  # full result dict
                    })

    fine_results.sort(key=lambda x: x['savings'], reverse=True)

    if save and fine_results:
        path = os.path.join(DATA_DIR, 'scan_2026_fine.json')
        # Don't serialize full 'result' to JSON (too large, has numpy arrays)
        save_results = []
        for r in fine_results:
            sr = {k: v for k, v in r.items() if k != 'result'}
            save_results.append(sr)
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(save_results, f, indent=2, default=str)
        print(f"  Saved: {path}")

    return fine_results


# ============================================================
#  Phase 3: Local refinement (golden-section search on rp)
# ============================================================

def refine_best(fine_results, iterations=10, save=True):
    """
    Refine the single best solution using localized search.

    Uses a simple Nelder-Mead style search on (day, rm, rp).
    Since day is discrete and rm/rp are continuous in ranges,
    this searches (rm, rp) for the best day and its neighbors.

    Returns
    -------
    dict: Optimal trajectory solution
    """
    if not fine_results:
        print("[Refine] No fine results to refine.")
        return None

    best = fine_results[0]
    print(f"[Refine] Starting from: day={best['day']}, "
          f"rm={best['rm']:.0f}, rp={best['rp_target_AU']:.4f} AU, "
          f"dv={best['dv_total']:.3f} km/s")

    best_day = best['day']
    best_rm = best['rm']
    best_rp = best['rp_target']
    best_dv = best['dv_total']

    # Search neighbor days
    day_candidates = list(range(max(0, best_day - 3), min(365, best_day + 4)))

    # Local (rm, rp) grid around best
    for day in day_candidates:
        for i in range(iterations):
            # Perturb rm and rp
            rm_scale = 1.0 + 0.3 * (np.random.random() - 0.5)  # +/- 15%
            rp_scale = 1.0 + 0.2 * (np.random.random() - 0.5)  # +/- 10%

            rm_try = np.clip(best_rm * rm_scale, C_RM_MIN, C_RM_MAX)
            rp_try = np.clip(best_rp * rp_scale, C_RP_MIN, C_RP_MAX)

            best_sol = solve_best_trajectory(int(day), rm_try, rp_try)
            if best_sol is None:
                continue

            direct = direct_launch_delta_v(rp_try)
            savings = direct['delta_v_total'] - best_sol['delta_v_total']

            if savings > 0 and best_sol['delta_v_total'] < best_dv:
                best_dv = best_sol['delta_v_total']
                best_day = int(day)
                best_rm = rm_try
                best_rp = rp_try
                print(f"    Improved: day={best_day}, rm={best_rm:.0f}, "
                      f"rp={best_rp/AU:.4f} AU, dv={best_dv:.3f} km/s")

    # Final solve with best params
    final = solve_best_trajectory(int(best_day), best_rm, best_rp)
    if final is None:
        print("[Refine] Final solve failed!")
        return None

    direct = direct_launch_delta_v(best_rp)
    optimal = {
        't0_day': int(best_day),
        'rm_optimal': best_rm,
        'rp_optimal': best_rp,
        'rp_optimal_AU': best_rp / AU,
        'delta_v_launch': final['delta_v_launch'],
        'delta_v_reentry': final['delta_v_reentry'],
        'delta_v_total': final['delta_v_total'],
        'savings': direct['delta_v_total'] - final['delta_v_total'],
        'savings_pct': 100.0 * (direct['delta_v_total'] - final['delta_v_total']) / direct['delta_v_total'],
        'delta_v_direct': direct['delta_v_total'],
        'side': final.get('side', 'N/A'),
        'trajectory': final,
    }

    print(f"\n[Optimal Solution]")
    print(f"  t0*: day {optimal['t0_day']} of 2026")
    print(f"  rm*: {optimal['rm_optimal']:.0f} km")
    print(f"  rp*: {optimal['rp_optimal_AU']:.4f} AU")
    print(f"  dv_total: {optimal['delta_v_total']:.3f} km/s")
    print(f"  Savings vs direct: {optimal['savings']:.3f} km/s ({optimal['savings_pct']:.2f}%)")
    print(f"  Side: {optimal['side']}")

    if save:
        path = os.path.join(DATA_DIR, 'scan_2026_optimal.json')
        save_dict = {k: v for k, v in optimal.items() if k != 'trajectory'}
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(save_dict, f, indent=2, default=str)
        print(f"  Saved: {path}")

    return optimal


# ============================================================
#  Full pipeline
# ============================================================

def find_optimal_2026(mode='coarse'):
    """
    Full optimization pipeline.

    Parameters
    ----------
    mode : str
        'quick'  - 28-day demo scan
        'coarse' - Phase 1 only (fastest)
        'full'   - Phase 1 + 2 + 3 (complete)

    Returns
    -------
    dict: Final optimal solution
    """
    print("=" * 60)
    print(f"  M6: Optimal Launch Window Scanner (mode={mode})")
    print("=" * 60)

    if mode == 'quick':
        n_days = 28
    else:
        n_days = 365

    # Phase 1: Coarse
    coarse = coarse_scan(n_days=n_days, save=True)

    if mode == 'coarse' or mode == 'quick':
        print("\n[M6] Coarse scan complete. Use --fine for full optimization.")
        return coarse.get('global_best')

    # Phase 2: Fine
    fine = fine_scan(coarse, top_n=10, save=True)

    # Phase 3: Refine
    optimal = refine_best(fine, save=True)

    print("\n[M6] Full optimization complete!")
    return optimal


# ============================================================
#  Utility: Generate dv_vs_day curve for plotting
# ============================================================

def generate_dv_curve(rp=0.2 * AU, n_days=365, save=True):
    """
    Generate Dv_total vs launch day curve for a fixed rp.
    Used by M8 (visualization).
    """
    print(f"[dv_curve] Scanning {n_days} days for rp={rp/AU:.2f} AU...")
    days = []
    dv_assist = []
    dv_direct_val = []
    sides = []

    for day in range(n_days):
        best = solve_best_trajectory(day, 5000.0, rp)
        if best is not None:
            days.append(day)
            dv_assist.append(best['delta_v_total'])
            dv_direct_val.append(direct_launch_delta_v(rp)['delta_v_total'])
            sides.append(best.get('side', 'N/A'))

    data = {
        'rp_AU': rp / AU,
        'days': days,
        'dv_assist': dv_assist,
        'dv_direct': dv_direct_val if dv_direct_val else [direct_launch_delta_v(rp)['delta_v_total']],
        'sides': sides,
    }

    if save:
        path = os.path.join(DATA_DIR, 'dv_vs_day.json')
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"  Saved: {path}")

    return data


# ============================================================
#  Main
# ============================================================

if __name__ == '__main__':
    mode = 'coarse'
    if '--fine' in sys.argv:
        mode = 'full'
    elif '--quick' in sys.argv:
        mode = 'quick'

    result = find_optimal_2026(mode=mode)

    # Also generate dv curve for rp=0.2 AU
    if mode in ('coarse', 'quick'):
        print("\n--- Dv vs Day curve ---")
        generate_dv_curve(rp=0.2 * AU, n_days=28 if mode == 'quick' else 365)
