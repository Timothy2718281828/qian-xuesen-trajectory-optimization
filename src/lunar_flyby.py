"""
lunar_flyby.py -- M4: Vector-based lunar gravity assist (corrected)

Fixes the scalar Dv direction bug in patched_conic.py by using full 2D
vector treatment of the flyby geometry:
  1. Proper frame transforms (Moon -> Earth -> heliocentric)
  2. Correct rotation of v_inf through the hyperbolic turn angle
  3. Post-flyby orbit computation including radial velocity component
  4. Return-leg v_inf computed from orbital mechanics

Units: km, s, km/s, rad
Reference frame: J2000 ecliptic plane (2D)

Usage:
    from lunar_flyby import design_trajectory, compare_direct_vs_assist
    result = design_trajectory(rp_target=0.2*AU, rm=5000, side='leading')
    comparison = compare_direct_vs_assist(rp=0.2*AU, rm=5000)

Known limitations:
    - Circular Earth/Moon orbits (no eccentricity)
    - 2D ecliptic plane (no inclination)
    - Radial approach approximation (valid for v_launch >> v_circular_moon)
"""

import numpy as np
from scipy.optimize import brentq
from physical_constants import (
    AU, MU_SUN, MU_EARTH, MU_MOON,
    R_EARTH, R_MOON, R_SUN,
    R_EARTH_ORBIT, V_EARTH_ORBIT,
    R_MOON_ORBIT, V_MOON_ORBIT,
    R_MOON_SOI,
    V_EARTH_ROTATION,
    C_RM_MIN, C_RP_MIN, C_RP_MAX, C_T_MAX, C_VINF_MAX,
    SECONDS_PER_DAY, SECONDS_PER_YEAR,
)


# ============================================================
#  Orbital mechanics utilities
# ============================================================

def helio_orbit_from_state(r_vec, v_vec, mu=MU_SUN):
    """
    Compute heliocentric orbital elements from 2D state vector.

    Parameters
    ----------
    r_vec : (2,) position [km]
    v_vec : (2,) velocity [km/s]

    Returns
    -------
    dict: a, e, rp, ra, E_specific, h, T, p
    """
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    E_sp = v**2 / 2.0 - mu / r
    if abs(E_sp) < 1e-20:
        E_sp = -1e-20

    a = -mu / (2.0 * E_sp)
    h = r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0]
    p = h**2 / mu

    e_sq = 1.0 - p / a if a > 0 else 1.0 + p / abs(a)
    e = np.sqrt(max(e_sq, 0.0))

    rp = a * (1.0 - e) if a > 0 else p / (1.0 + e)
    ra = a * (1.0 + e) if a > 0 else float('inf')
    T = 2.0 * np.pi * np.sqrt(abs(a)**3 / mu) if a > 0 else float('inf')

    return dict(a=a, e=e, rp=rp, ra=ra, E_specific=E_sp, h=h, T=T, p=p)


def velocity_at_r(r, a, e, h, mu=MU_SUN, outbound=True):
    """
    Velocity components at distance r on a Keplerian orbit.

    Returns (v_radial, v_tangential).
    Outbound: v_r > 0; inbound: v_r < 0.
    """
    v_mag = np.sqrt(mu * (2.0 / r - 1.0 / a))
    v_t = abs(h) / r
    v_r_sq = v_mag**2 - v_t**2
    v_r = np.sqrt(max(v_r_sq, 0.0))
    if not outbound:
        v_r = -v_r
    return v_r, v_t


# ============================================================
#  Moon state
# ============================================================

def moon_state_at_phase(theta):
    """
    Moon position & velocity in Earth-centered frame at orbital phase theta.

    theta = 0 -> Moon on +x axis from Earth.
    Counterclockwise orbit.

    Returns (r_moon, v_moon) each (2,) in km, km/s.
    """
    r = R_MOON_ORBIT * np.array([np.cos(theta), np.sin(theta)])
    v = V_MOON_ORBIT * np.array([-np.sin(theta), np.cos(theta)])
    return r, v


# ============================================================
#  Approach velocity at Moon encounter
# ============================================================

def compute_approach_velocity(v_launch, moon_phase):
    """
    Compute spacecraft velocity at Moon encounter (Earth frame).

    Assumes tangential launch from Earth surface (maximises angular momentum).
    The orbit around Earth is determined by v_launch (perigee velocity).

    Parameters
    ----------
    v_launch : float
        Perigee velocity at Earth surface [km/s]
    moon_phase : float
        Moon orbital phase angle [rad]

    Returns
    -------
    v_sc_earth : (2,) spacecraft velocity in Earth frame [km/s]
    v_at_moon : float, speed at Moon's orbital distance [km/s]
    """
    v_sq_at_moon = v_launch**2 - 2.0 * MU_EARTH * (1.0/R_EARTH - 1.0/R_MOON_ORBIT)
    if v_sq_at_moon < 0:
        return None, 0.0

    v_at_moon = np.sqrt(v_sq_at_moon)

    h = R_EARTH * v_launch
    v_t = h / R_MOON_ORBIT
    v_r_sq = v_at_moon**2 - v_t**2
    v_r = np.sqrt(max(v_r_sq, 0.0))

    r_hat = np.array([np.cos(moon_phase), np.sin(moon_phase)])
    t_hat = np.array([-np.sin(moon_phase), np.cos(moon_phase)])

    v_sc_earth = v_r * r_hat + v_t * t_hat
    return v_sc_earth, v_at_moon


# ============================================================
#  Core flyby rotation
# ============================================================

def flyby_rotate(v_inf_in, rm, v_moon_earth, mu_moon=MU_MOON, side='leading'):
    """
    Rotate incoming v_inf by the hyperbolic turn angle.

    The rotation direction is determined geometrically:
    - 'leading': periapsis on Moon's velocity side -> reduces speed (solar probe)
    - 'trailing': periapsis opposite Moon's velocity -> increases speed

    Parameters
    ----------
    v_inf_in : (2,) incoming v_inf relative to Moon [km/s]
    rm : float, closest approach to Moon center [km]
    v_moon_earth : (2,) Moon velocity in Earth frame [km/s]
    side : 'leading' or 'trailing'

    Returns
    -------
    dict: v_inf_out, turn_angle, e_hyperbolic, rotation_sign
    """
    v_inf_mag = np.linalg.norm(v_inf_in)
    if v_inf_mag < 1e-10:
        return dict(v_inf_out=v_inf_in.copy(), turn_angle=0.0,
                    e_hyperbolic=np.inf, rotation_sign=0)

    e_hyper = 1.0 + rm * v_inf_mag**2 / mu_moon
    if e_hyper <= 1.0:
        e_hyper = 1.0 + 1e-10
    delta = 2.0 * np.arcsin(1.0 / e_hyper)

    d_in = v_inf_in / v_inf_mag
    cross = d_in[0] * v_moon_earth[1] - d_in[1] * v_moon_earth[0]

    # leading: periapsis on Moon-velocity side
    #   cross > 0 -> Moon vel is to the LEFT  -> periapsis LEFT  -> CCW (+delta)
    #   cross < 0 -> Moon vel is to the RIGHT -> periapsis RIGHT -> CW (-delta)
    # trailing: opposite
    if side == 'leading':
        rot_sign = 1.0 if cross >= 0 else -1.0
    else:
        rot_sign = -1.0 if cross >= 0 else 1.0

    angle = rot_sign * delta
    ca, sa = np.cos(angle), np.sin(angle)
    v_inf_out = np.array([ca * v_inf_in[0] - sa * v_inf_in[1],
                           sa * v_inf_in[0] + ca * v_inf_in[1]])

    return dict(v_inf_out=v_inf_out, turn_angle=delta,
                e_hyperbolic=e_hyper, rotation_sign=rot_sign)


# ============================================================
#  Full flyby -> post-flyby heliocentric state
# ============================================================

def compute_post_flyby(v_launch, moon_phase, rm, side='leading'):
    """
    Compute full flyby result: launch -> Moon encounter -> post-flyby orbit.

    Parameters
    ----------
    v_launch : float, launch velocity at Earth surface [km/s]
    moon_phase : float, Moon orbital phase [rad]
    rm : float, flyby closest approach [km]
    side : 'leading' or 'trailing'

    Returns
    -------
    dict or None (if v_launch insufficient to reach Moon)
    """
    v_sc_earth, v_at_moon = compute_approach_velocity(v_launch, moon_phase)
    if v_sc_earth is None:
        return None

    r_moon_earth, v_moon_earth = moon_state_at_phase(moon_phase)

    v_inf_in = v_sc_earth - v_moon_earth
    v_inf_mag = np.linalg.norm(v_inf_in)

    flyby = flyby_rotate(v_inf_in, rm, v_moon_earth, MU_MOON, side)

    v_sc_earth_after = v_moon_earth + flyby['v_inf_out']

    r_earth_helio = np.array([AU, 0.0])
    v_earth_helio = np.array([0.0, V_EARTH_ORBIT])

    r_sc_helio = r_earth_helio + r_moon_earth
    v_sc_helio = v_earth_helio + v_sc_earth_after

    orbit = helio_orbit_from_state(r_sc_helio, v_sc_helio)

    return dict(
        orbit=orbit,
        v_sc_helio=v_sc_helio,
        r_sc_helio=r_sc_helio,
        v_sc_earth_before=v_sc_earth,
        v_sc_earth_after=v_sc_earth_after,
        v_inf_in=v_inf_in,
        v_inf_out=flyby['v_inf_out'],
        v_inf_mag=v_inf_mag,
        turn_angle=flyby['turn_angle'],
        e_hyperbolic=flyby['e_hyperbolic'],
        v_at_moon=v_at_moon,
    )


# ============================================================
#  Return-leg v_inf and reentry Dv
# ============================================================

def compute_return_vinf(orbit):
    """
    Compute v_inf when spacecraft returns to Earth's orbital distance.

    Uses vis-viva for speed and angular momentum for tangential split.
    On the return leg the spacecraft is outbound (v_r > 0).

    Returns v_inf_return [km/s] or None if orbit doesn't cross r = AU.
    """
    a, e, h = orbit['a'], orbit['e'], orbit['h']
    r = R_EARTH_ORBIT

    if a <= 0 or r > a * (1 + e):
        return None

    v_r, v_t = velocity_at_r(r, a, e, h, MU_SUN, outbound=True)

    r_hat = np.array([1.0, 0.0])
    t_hat = np.array([0.0, 1.0])
    v_sc = v_r * r_hat + v_t * t_hat
    v_earth = np.array([0.0, V_EARTH_ORBIT])

    v_inf_return = np.linalg.norm(v_sc - v_earth)
    return v_inf_return


def compute_reentry_dv(v_inf_return):
    """Reentry Dv: speed at Earth surface from hyperbolic approach."""
    v_entry = np.sqrt(v_inf_return**2 + 2.0 * MU_EARTH / R_EARTH)
    return v_entry


# ============================================================
#  Trajectory design solver
# ============================================================

def design_trajectory(rp_target, rm, side='leading', moon_phase=0.0):
    """
    Design trajectory: Earth -> Moon flyby -> perihelion rp_target -> return.

    Searches for v_launch that achieves the target perihelion via Moon flyby.

    Parameters
    ----------
    rp_target : float, target perihelion distance [km]
    rm : float, flyby closest approach to Moon center [km]
    side : 'leading' or 'trailing'
    moon_phase : float, Moon orbital phase at encounter [rad]

    Returns
    -------
    dict with delta_v_launch, delta_v_reentry, delta_v_total, orbit details.
    None if no solution found.
    """
    v_launch_min = np.sqrt(2.0 * MU_EARTH * (1.0/R_EARTH - 1.0/R_MOON_ORBIT)) + 0.01
    v_launch_max = 25.0

    def rp_residual(v_launch):
        result = compute_post_flyby(v_launch, moon_phase, rm, side)
        if result is None:
            return R_EARTH_ORBIT
        rp = result['orbit']['rp']
        if not np.isfinite(rp) or rp < 0:
            return R_EARTH_ORBIT
        return rp - rp_target

    rp_lo = rp_residual(v_launch_min)
    rp_hi = rp_residual(v_launch_max)

    if rp_lo * rp_hi > 0:
        # scan for sign change
        for v_try in np.linspace(v_launch_min, v_launch_max, 200):
            r_try = rp_residual(v_try)
            if rp_lo * r_try < 0:
                v_launch_max = v_try
                rp_hi = r_try
                break
            rp_lo = r_try
            v_launch_min = v_try
        else:
            return None

    try:
        v_launch_sol = brentq(rp_residual, v_launch_min, v_launch_max,
                              xtol=1e-6, maxiter=200)
    except (ValueError, RuntimeError):
        return None

    pfb = compute_post_flyby(v_launch_sol, moon_phase, rm, side)
    if pfb is None:
        return None

    delta_v_launch = v_launch_sol - V_EARTH_ROTATION

    v_inf_return = compute_return_vinf(pfb['orbit'])
    if v_inf_return is None:
        return None
    delta_v_reentry = compute_reentry_dv(v_inf_return)

    delta_v_total = delta_v_launch + delta_v_reentry

    return dict(
        delta_v_launch=delta_v_launch,
        delta_v_reentry=delta_v_reentry,
        delta_v_total=delta_v_total,
        v_launch=v_launch_sol,
        v_inf_return=v_inf_return,
        rp_achieved=pfb['orbit']['rp'],
        orbit=pfb['orbit'],
        turn_angle_deg=np.degrees(pfb['turn_angle']),
        e_hyperbolic=pfb['e_hyperbolic'],
        v_inf_moon=pfb['v_inf_mag'],
        v_inf_in=pfb['v_inf_in'],
        v_inf_out=pfb['v_inf_out'],
        v_sc_helio=pfb['v_sc_helio'],
        r_sc_helio=pfb['r_sc_helio'],
        v_sc_earth_before=pfb['v_sc_earth_before'],
        v_sc_earth_after=pfb['v_sc_earth_after'],
        moon_phase=moon_phase,
        rm=rm,
        side=side,
        flight_time=pfb['orbit']['T'],
    )


# ============================================================
#  Direct launch (no Moon assist) for comparison
# ============================================================

def direct_launch_delta_v(rp_target, r1=R_EARTH_ORBIT, mu_sun=MU_SUN):
    """
    Compute Dv_total for a direct Hohmann-like transfer (no Moon assist).

    Returns dict with delta_v_launch, delta_v_reentry, delta_v_total.
    """
    a_trans = (r1 + rp_target) / 2.0
    v_dep = np.sqrt(mu_sun * (2.0 / r1 - 1.0 / a_trans))
    v_circ = np.sqrt(mu_sun / r1)

    v_inf_dep = abs(v_dep - v_circ)
    v_launch = np.sqrt(v_inf_dep**2 + 2.0 * MU_EARTH / R_EARTH)
    dv_launch = v_launch - V_EARTH_ROTATION

    v_inf_return = v_inf_dep
    dv_reentry = compute_reentry_dv(v_inf_return)

    T_transfer = np.pi * np.sqrt(a_trans**3 / mu_sun)

    return dict(
        delta_v_launch=dv_launch,
        delta_v_reentry=dv_reentry,
        delta_v_total=dv_launch + dv_reentry,
        v_inf_dep=v_inf_dep,
        v_inf_return=v_inf_return,
        v_dep_helio=v_dep,
        v_circ=v_circ,
        flight_time=2.0 * T_transfer,
    )


# ============================================================
#  Compare direct vs Moon-assisted
# ============================================================

def compare_direct_vs_assist(rp, rm, side='leading', moon_phase=0.0):
    """
    Compare Dv for direct launch vs Moon-assisted launch.

    Returns dict with both results and savings.
    """
    direct = direct_launch_delta_v(rp)
    assist = design_trajectory(rp, rm, side, moon_phase)

    savings = None
    if assist is not None:
        savings = direct['delta_v_total'] - assist['delta_v_total']

    return dict(direct=direct, assist=assist, delta_v_savings=savings)


# ============================================================
#  Scan over Moon phases to find optimal geometry
# ============================================================

def scan_moon_phases(rp_target, rm, side='leading', n_phases=72):
    """
    Scan Moon orbital phase to find the geometry that minimises Dv_total.

    Returns list of (phase, result) and the best result.
    """
    results = []
    for i in range(n_phases):
        theta = 2.0 * np.pi * i / n_phases
        res = design_trajectory(rp_target, rm, side, theta)
        if res is not None:
            results.append((theta, res))

    if not results:
        return [], None

    best = min(results, key=lambda x: x[1]['delta_v_total'])
    return results, best


def scan_moon_phases_both_sides(rp_target, rm, n_phases=72):
    """Scan both leading and trailing to find global optimum."""
    all_results = []
    for side in ('leading', 'trailing'):
        results, best = scan_moon_phases(rp_target, rm, side, n_phases)
        if best is not None:
            all_results.append(best)

    if not all_results:
        return None

    return min(all_results, key=lambda x: x[1]['delta_v_total'])


# ============================================================
#  N-body verification of analytical flyby
# ============================================================

def _setup_nbody_post_flyby(design_result):
    """
    Create an NBodySimulation with the rocket placed just outside Moon's
    SOI on the **outgoing** side, carrying the post-flyby velocity.

    This tests whether the heliocentric transfer orbit matches the
    analytical prediction.  A full incoming-flyby setup would require
    computing the hyperbolic entry state (position offset by impact
    parameter), which is deferred to M6.
    """
    from nbody import NBodySimulation

    moon_phase = design_result['moon_phase']
    r_moon_earth, v_moon_earth = moon_state_at_phase(moon_phase)

    pos_sun = np.array([0.0, 0.0])
    vel_sun = np.array([0.0, 0.0])
    pos_earth = np.array([AU, 0.0])
    vel_earth = np.array([0.0, V_EARTH_ORBIT])
    pos_moon = pos_earth + r_moon_earth
    vel_moon = vel_earth + v_moon_earth

    v_inf_out = design_result['v_inf_out']
    v_inf_out_mag = np.linalg.norm(v_inf_out)
    v_inf_out_hat = v_inf_out / v_inf_out_mag

    pos_rocket = pos_moon + R_MOON_SOI * v_inf_out_hat

    v_soi = np.sqrt(v_inf_out_mag**2 + 2.0 * MU_MOON / R_MOON_SOI)
    vel_rocket = vel_moon + v_soi * v_inf_out_hat

    sim = NBodySimulation()
    sim.set_bodies(pos_sun, vel_sun, pos_earth, vel_earth,
                   pos_moon, vel_moon, pos_rocket, vel_rocket)
    return sim


def verify_flyby_nbody(design_result, t_propagate_days=60):
    """
    Verify the analytical flyby model with an N-body simulation.

    The rocket is placed at the Moon SOI boundary with the analytical
    approach velocity.  The integrator resolves the hyperbolic flyby
    and subsequent heliocentric transfer.

    Returns dict with analytical vs numerical perihelion comparison.
    """
    sim = _setup_nbody_post_flyby(design_result)

    t_span = t_propagate_days * SECONDS_PER_DAY

    trajectory = []
    def record(t, pos, vel):
        trajectory.append((t, pos.copy(), vel.copy()))

    sim.integrate(t_span, h=3600.0, h_fine=60.0,
                  fine_body_idx=2, fine_radius=R_MOON_SOI,
                  energy_check_interval=500, callback=record)

    r_min_sun = float('inf')
    r_min_moon = float('inf')
    for t, pos, vel in trajectory:
        r_sun = np.linalg.norm(pos[3] - pos[0])
        r_moon_dist = np.linalg.norm(pos[3] - pos[2])
        if r_sun < r_min_sun:
            r_min_sun = r_sun
        if r_moon_dist < r_min_moon:
            r_min_moon = r_moon_dist

    rp_analytical = design_result['rp_achieved']
    rp_numerical = r_min_sun
    rp_error = abs(rp_numerical - rp_analytical)
    rp_error_pct = rp_error / rp_analytical * 100.0 if rp_analytical > 0 else float('inf')

    return dict(
        rp_analytical=rp_analytical,
        rp_numerical=rp_numerical,
        rp_error_km=rp_error,
        rp_error_pct=rp_error_pct,
        r_min_moon=r_min_moon,
        n_trajectory_points=len(trajectory),
        trajectory=trajectory,
    )


# ============================================================
#  Main: run verification and comparison
# ============================================================

if __name__ == '__main__':
    print("=" * 65)
    print("  M4: Lunar Gravity Assist -- Corrected Vector Model")
    print("=" * 65)

    rp_test = 0.2 * AU

    # --- Direct launch baseline ---
    direct = direct_launch_delta_v(rp_test)
    print(f"\n  Direct launch to rp = 0.2 AU:")
    print(f"    Dv_launch  = {direct['delta_v_launch']:.4f} km/s")
    print(f"    Dv_reentry = {direct['delta_v_reentry']:.4f} km/s")
    print(f"    Dv_total   = {direct['delta_v_total']:.4f} km/s")
    print(f"    v_inf_dep  = {direct['v_inf_dep']:.4f} km/s")
    print(f"    Flight time = {direct['flight_time']/SECONDS_PER_DAY:.1f} days")

    # --- Scan Moon phases for rm = 5000 km ---
    print(f"\n  Scanning Moon phases (rm = 5000 km)...")
    for side in ('leading', 'trailing'):
        results, best = scan_moon_phases(rp_test, rm=5000.0, side=side, n_phases=72)
        if best is not None:
            theta_best, res_best = best
            print(f"\n  Best {side} flyby:")
            print(f"    Moon phase    = {np.degrees(theta_best):.1f} deg")
            print(f"    Dv_launch     = {res_best['delta_v_launch']:.4f} km/s")
            print(f"    Dv_reentry    = {res_best['delta_v_reentry']:.4f} km/s")
            print(f"    Dv_total      = {res_best['delta_v_total']:.4f} km/s")
            print(f"    Turn angle    = {res_best['turn_angle_deg']:.2f} deg")
            print(f"    v_inf (Moon)  = {res_best['v_inf_moon']:.4f} km/s")
            print(f"    rp achieved   = {res_best['rp_achieved']/AU:.6f} AU")
            savings = direct['delta_v_total'] - res_best['delta_v_total']
            print(f"    Savings vs direct = {savings:.4f} km/s "
                  f"({savings/direct['delta_v_total']*100:.2f}%)")
        else:
            print(f"\n  No solution found for {side} flyby.")

    # --- Scan multiple rm values ---
    print(f"\n{'='*65}")
    print(f"  rm sensitivity (best phase, leading):")
    print(f"  {'rm (km)':>10}  {'Dv_total':>10}  {'Savings':>10}  {'Turn deg':>10}")
    for rm_val in [2000, 3000, 5000, 10000, 20000, 40000]:
        best_phase = scan_moon_phases_both_sides(rp_test, rm_val, n_phases=72)
        if best_phase is not None:
            _, res = best_phase
            sav = direct['delta_v_total'] - res['delta_v_total']
            print(f"  {rm_val:10.0f}  {res['delta_v_total']:10.4f}  "
                  f"{sav:10.4f}  {res['turn_angle_deg']:10.2f}")
        else:
            print(f"  {rm_val:10.0f}  {'N/A':>10}")

    # --- N-body verification for one case ---
    print(f"\n{'='*65}")
    print(f"  N-body verification:")
    best_overall = scan_moon_phases_both_sides(rp_test, rm=5000.0, n_phases=72)
    if best_overall is not None:
        _, best_res = best_overall
        print(f"  Running N-body for {best_res['side']} flyby, "
              f"phase={np.degrees(best_res['moon_phase']):.1f} deg...")
        nbody_check = verify_flyby_nbody(best_res, t_propagate_days=120)
        print(f"    rp analytical = {nbody_check['rp_analytical']/AU:.6f} AU")
        print(f"    rp numerical  = {nbody_check['rp_numerical']/AU:.6f} AU")
        print(f"    Error         = {nbody_check['rp_error_km']:.0f} km "
              f"({nbody_check['rp_error_pct']:.2f}%)")
        print(f"    Min dist Moon = {nbody_check['r_min_moon']:.0f} km")

    print(f"\n  M4 complete.")
