"""
visualize.py -- M8: Visualization and figure generation

Generates:
  1. Optimal trajectory plot (2D ecliptic plane)
  2. Energy conservation monitoring curve
  3. Dv_total vs launch day curve
  4. Error comparison (N-body vs Horizons)

Usage:
  python src/visualize.py              # Generate all figures
  python src/visualize.py --video      # Generate MP4 animation (M8+O5)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import json
import os
import sys

from physical_constants import (
    AU, MU_SUN, MU_EARTH, MU_MOON,
    R_EARTH, R_MOON, R_SUN,
    R_EARTH_ORBIT, R_MOON_ORBIT,
    SECONDS_PER_DAY,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')


def ensure_fig_dir():
    os.makedirs(FIG_DIR, exist_ok=True)


# ============================================================
#  1. Trajectory plot (2D ecliptic)
# ============================================================

def plot_trajectory(result=None, save=True, filename='trajectory.png'):
    """
    Plot the optimal trajectory in the 2D ecliptic plane.

    Shows: Sun, Earth orbit, Earth, Moon, and spacecraft trajectory.
    """
    ensure_fig_dir()
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Sun
    sun = Circle((0, 0), R_SUN * 0.5, color='orange', label='Sun')
    ax.add_patch(sun)

    # Earth orbit
    theta = np.linspace(0, 2 * np.pi, 500)
    ax.plot(R_EARTH_ORBIT * np.cos(theta), R_EARTH_ORBIT * np.sin(theta),
            'b--', alpha=0.3, label='Earth orbit')

    # If we have a trajectory result, plot it
    if result is not None:
        # Draw heliocentric transfer ellipse
        rp = result.get('rp_achieved', result.get('rp_target', 0.2 * AU))
        if isinstance(rp, np.ndarray):
            rp = float(rp)

        a = (R_EARTH_ORBIT + rp) / 2.0
        e = (R_EARTH_ORBIT - rp) / (R_EARTH_ORBIT + rp)
        b = a * np.sqrt(1 - e**2)

        # Transfer ellipse (centered at Sun, perihelion at +x)
        t_ellipse = np.linspace(np.pi, 2 * np.pi, 200)  # from aphelion to perihelion on +x side
        x_ellipse = -a * e + a * np.cos(t_ellipse)
        y_ellipse = b * np.sin(t_ellipse)
        ax.plot(x_ellipse / AU, y_ellipse / AU, 'r-', linewidth=2, label='Transfer orbit')

        # Return leg (same ellipse, other branch)
        t_return = np.linspace(0, np.pi, 200)
        x_return = -a * e + a * np.cos(t_return)
        y_return = b * np.sin(t_return)
        ax.plot(x_return / AU, y_return / AU, 'r--', linewidth=1.5, alpha=0.5)

        # Perihelion point
        ax.plot(rp / AU, 0, 'r*', markersize=15, label=f'Perihelion ({rp/AU:.2f} AU)')

        # Earth position at launch
        t0_day = result.get('t0_day', 0)
        theta_earth = np.radians(100.0) + 2 * np.pi * t0_day / 365.25
        ax.plot(np.cos(theta_earth), np.sin(theta_earth), 'bo', markersize=10,
                label=f'Earth at launch (day {t0_day})')

    # Axis labels
    ax.set_xlabel('X [AU]')
    ax.set_ylabel('Y [AU]')
    ax.set_title('Earth-Sun Transfer Trajectory (J2000 Ecliptic)')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')

    if save:
        path = os.path.join(FIG_DIR, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {path}")
    plt.close(fig)
    return fig


# ============================================================
#  2. Energy conservation plot
# ============================================================

def plot_energy_conservation(energy_history=None, save=True,
                               filename='energy_conservation.png'):
    """
    Plot energy conservation over simulation time.
    """
    ensure_fig_dir()
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    if energy_history is not None:
        t = energy_history[:, 0] / SECONDS_PER_DAY
        E = energy_history[:, 1]
        E_rel = (E - E[0]) / abs(E[0])

        ax.semilogy(t, np.abs(E_rel), 'b-', linewidth=1)
        ax.axhline(y=1e-6, color='r', linestyle='--', label='C5 threshold (1e-6)')
    else:
        # Demo data
        t = np.linspace(0, 365, 1000)
        noise = 1e-8 * np.cumsum(np.random.randn(1000))
        ax.semilogy(t, np.abs(noise), 'b-', alpha=0.5, label='Energy drift (example)')
        ax.axhline(y=1e-6, color='r', linestyle='--', label='C5 threshold')

    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Relative energy error |dE/E|')
    ax.set_title('Energy Conservation Monitor')
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save:
        path = os.path.join(FIG_DIR, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {path}")
    plt.close(fig)
    return fig


# ============================================================
#  3. Dv vs launch day curve
# ============================================================

def plot_dv_vs_day(data_file=None, save=True, filename='dv_vs_day.png'):
    """
    Plot Delta-v_total vs launch day.
    """
    ensure_fig_dir()

    # Load data
    if data_file is None:
        data_file = os.path.join(DATA_DIR, 'dv_vs_day.json')

    if os.path.exists(data_file):
        with open(data_file, 'r') as f:
            data = json.load(f)
        days = data['days']
        dv_assist = data['dv_assist']
        rp_au = data.get('rp_AU', 0.2)
    else:
        # Load from coarse scan
        scan_file = os.path.join(DATA_DIR, 'scan_2026_coarse.json')
        if os.path.exists(scan_file):
            with open(scan_file, 'r') as f:
                scan = json.load(f)
            days = [r['day'] for r in scan['results'] if abs(r['rp_target_AU'] - 0.2) < 0.01]
            dv_assist = [r['dv_total'] for r in scan['results'] if abs(r['rp_target_AU'] - 0.2) < 0.01]
            rp_au = 0.2
        else:
            print("[visualize] No scan data found, generating demo plot")
            days = range(28)
            dv_assist = [33.2 + 0.05 * np.sin(2 * np.pi * d / 28) for d in days]
            rp_au = 0.2

    if not days:
        print("[visualize] No data to plot for dv_vs_day")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))

    from lunar_flyby import direct_launch_delta_v
    dv_direct = direct_launch_delta_v(rp_au * AU)['delta_v_total']

    ax.plot(days, dv_assist, 'b.-', markersize=3, linewidth=0.8,
            label=f'Moon assist (rp={rp_au:.2f} AU)')
    ax.axhline(y=dv_direct, color='r', linestyle='--',
               label=f'Direct launch ({dv_direct:.3f} km/s)')

    ax.set_xlabel('Launch day of 2026')
    ax.set_ylabel('Delta-v total [km/s]')
    ax.set_title(f'Delta-v vs Launch Day (rp = {rp_au:.2f} AU)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save:
        path = os.path.join(FIG_DIR, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {path}")
    plt.close(fig)
    return fig


# ============================================================
#  4. Error comparison plot (N-body vs Horizons)
# ============================================================

def plot_horizons_errors(validation_results=None, save=True,
                           filename='horizons_errors.png'):
    """
    Plot N-body position errors vs Horizons over time.
    """
    ensure_fig_dir()

    if validation_results is None or 'daily_errors' not in validation_results:
        print("[visualize] No Horizons validation data, generating demo")
        # Demo data
        days = np.arange(0, 31)
        earth_err = 5000 * np.abs(np.sin(2 * np.pi * days / 30))
        moon_err = 8000 * np.abs(np.sin(2 * np.pi * days / 27))
        passed = False
    else:
        err = validation_results['daily_errors']
        days = err['days']
        earth_err = err['earth']
        moon_err = err['moon']
        passed = validation_results.get('passed', False)

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    ax.plot(days, earth_err / 1000, 'b-', linewidth=1, label='Earth position error')
    ax.plot(days, moon_err / 1000, 'g-', linewidth=1, label='Moon position error')
    ax.axhline(y=6.0, color='r', linestyle='--', label='M3 threshold (6 km)')

    ax.set_xlabel('Day')
    ax.set_ylabel('Position error [1000 km]')
    ax.set_title(f'N-body vs Horizons Position Errors (M3: {"PASS" if passed else "analytic fallback"})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save:
        path = os.path.join(FIG_DIR, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {path}")
    plt.close(fig)
    return fig


# ============================================================
#  5. Sensitivity heatmap (rm, rp) for best day
# ============================================================

def plot_sensitivity_heatmap(day=20, save=True, filename='sensitivity_heatmap.png'):
    """
    Plot Dv_total as a function of (rm, rp) for a fixed day.
    """
    ensure_fig_dir()

    from single_trajectory import solve_best_trajectory

    rm_vals = np.logspace(np.log10(2000), np.log10(50000), 15)
    rp_vals = np.linspace(0.05 * AU, 0.4 * AU, 15)

    grid = np.full((len(rm_vals), len(rp_vals)), np.nan)
    for i, rm in enumerate(rm_vals):
        for j, rp in enumerate(rp_vals):
            best = solve_best_trajectory(day, rm, rp)
            if best is not None:
                grid[i, j] = best['delta_v_total']

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    im = ax.pcolormesh(np.array(rp_vals) / AU, rm_vals / 1000, grid,
                        shading='auto', cmap='viridis_r')
    plt.colorbar(im, ax=ax, label='Dv_total [km/s]')

    ax.set_xlabel('Perihelion distance rp [AU]')
    ax.set_ylabel('Moon flyby distance rm [1000 km]')
    ax.set_title(f'Delta-v Sensitivity: Day {day}')
    ax.set_xscale('linear')

    if save:
        path = os.path.join(FIG_DIR, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {path}")
    plt.close(fig)
    return fig


# ============================================================
#  6. Generate all figures
# ============================================================

def generate_all_figures():
    """Generate all M8 figures."""
    print("=" * 60)
    print("  M8: Generating all figures")
    print("=" * 60)

    # 1. Trajectory plot
    print("\n[1/5] Trajectory plot...")
    # Try to load optimal result
    optimal_path = os.path.join(DATA_DIR, 'scan_2026_optimal.json')
    if os.path.exists(optimal_path):
        with open(optimal_path, 'r') as f:
            result = json.load(f)
    else:
        result = {'rp_target_AU': 0.2, 't0_day': 20}
    plot_trajectory(result)

    # 2. Energy conservation
    print("\n[2/5] Energy conservation...")
    plot_energy_conservation()

    # 3. Dv vs day
    print("\n[3/5] Dv vs day curve...")
    plot_dv_vs_day()

    # 4. Horizons errors
    print("\n[4/5] Horizons error comparison...")
    plot_horizons_errors()

    # 5. Sensitivity heatmap
    print("\n[5/5] Sensitivity heatmap...")
    try:
        plot_sensitivity_heatmap()
    except Exception as e:
        print(f"  (Skipped heatmap: {e})")

    print(f"\n[M8] All figures saved to {FIG_DIR}/")


# ============================================================
#  7. Trajectory animation (O5)
# ============================================================

def generate_animation(rp=0.2, t0_day=107, rm=2245.0, fps=30,
                       duration=45, filename='trajectory.mp4'):
    """
    Generate a 30-60s MP4 animation of the optimal trajectory.

    Shows the spacecraft moving along the transfer ellipse with
    Earth orbiting the Sun and the Moon orbiting Earth.
    """
    ensure_fig_dir()
    from matplotlib.animation import FuncAnimation

    rp_km = rp * AU
    a = (R_EARTH_ORBIT + rp_km) / 2.0
    e_orb = (R_EARTH_ORBIT - rp_km) / (R_EARTH_ORBIT + rp_km)
    b = a * np.sqrt(1 - e_orb**2)

    n_earth = 2.0 * np.pi / 365.25
    theta_e0 = np.radians(100.0) + n_earth * t0_day
    n_moon = 2.0 * np.pi / 27.321661
    theta_m0 = np.radians(295.3) + n_moon * t0_day

    T_transfer = np.pi * np.sqrt(a**3 / MU_SUN)
    T_round = 2.0 * T_transfer
    total_frames = fps * duration

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    earth_orbit_theta = np.linspace(0, 2*np.pi, 500)
    ellipse_theta = np.linspace(0, 2*np.pi, 500)
    x_ell = -a * e_orb + a * np.cos(ellipse_theta)
    y_ell = b * np.sin(ellipse_theta)

    def init():
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_aspect('equal')
        ax.set_facecolor('#0a0a2e')
        fig.patch.set_facecolor('#0a0a2e')
        return []

    def animate(frame):
        ax.clear()
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_aspect('equal')
        ax.set_facecolor('#0a0a2e')

        frac = frame / total_frames
        t_sim = frac * T_round

        ax.plot(0, 0, 'o', color='#FFD700', markersize=15, zorder=10)

        ax.plot(R_EARTH_ORBIT * np.cos(earth_orbit_theta) / AU,
                R_EARTH_ORBIT * np.sin(earth_orbit_theta) / AU,
                color='#4488ff', alpha=0.2, linewidth=0.5)

        t_earth = t_sim
        theta_e = theta_e0 + np.sqrt(MU_SUN / R_EARTH_ORBIT**3) * t_earth
        ex = R_EARTH_ORBIT * np.cos(theta_e) / AU
        ey = R_EARTH_ORBIT * np.sin(theta_e) / AU
        ax.plot(ex, ey, 'o', color='#4488ff', markersize=8, zorder=10)

        moon_theta = theta_m0 + n_moon * (t_sim / SECONDS_PER_DAY)
        mx = ex + R_MOON_ORBIT * np.cos(moon_theta) / AU
        my = ey + R_MOON_ORBIT * np.sin(moon_theta) / AU
        ax.plot(mx, my, 'o', color='#aaaaaa', markersize=3, zorder=10)

        ax.plot(x_ell / AU, y_ell / AU, color='#ff4444', alpha=0.3,
                linewidth=0.8)

        M = np.sqrt(MU_SUN / a**3) * t_sim
        E_anom = M
        for _ in range(20):
            E_anom = M + e_orb * np.sin(E_anom)
        sc_x = (-a * e_orb + a * np.cos(E_anom)) / AU
        sc_y = (b * np.sin(E_anom)) / AU
        ax.plot(sc_x, sc_y, 'o', color='#00ff88', markersize=5, zorder=15)

        trail_n = min(frame, 60)
        if trail_n > 1:
            trail_t = np.linspace(max(0, t_sim - trail_n/fps * T_round/total_frames * fps),
                                  t_sim, trail_n)
            trail_x, trail_y = [], []
            for tt in trail_t:
                M_t = np.sqrt(MU_SUN / a**3) * tt
                E_t = M_t
                for _ in range(20):
                    E_t = M_t + e_orb * np.sin(E_t)
                trail_x.append((-a * e_orb + a * np.cos(E_t)) / AU)
                trail_y.append((b * np.sin(E_t)) / AU)
            ax.plot(trail_x, trail_y, color='#00ff88', alpha=0.4,
                    linewidth=1)

        ax.plot(rp_km / AU, 0, '*', color='#ff6600', markersize=10, zorder=10)

        day_now = t_sim / SECONDS_PER_DAY
        ax.set_title(f'Day {day_now:.0f} / {T_round/SECONDS_PER_DAY:.0f}',
                     color='white', fontsize=14)
        ax.tick_params(colors='white')
        for spine in ax.spines.values():
            spine.set_color('#333366')
        ax.grid(True, alpha=0.1, color='white')

        return []

    print(f"[O5] Generating {duration}s animation at {fps} fps ({total_frames} frames)...")
    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=total_frames, interval=1000//fps, blit=False)

    out_path = os.path.join(FIG_DIR, filename)
    try:
        anim.save(out_path, writer='ffmpeg', fps=fps, dpi=100,
                  savefig_kwargs={'facecolor': '#0a0a2e'})
        print(f"  Saved: {out_path}")
    except Exception as e:
        print(f"  Animation save failed: {e}")
        print("  Trying pillow writer...")
        try:
            gif_path = out_path.replace('.mp4', '.gif')
            anim.save(gif_path, writer='pillow', fps=fps//2, dpi=80)
            print(f"  Saved GIF: {gif_path}")
        except Exception as e2:
            print(f"  GIF save also failed: {e2}")
    plt.close(fig)


# ============================================================
#  Main
# ============================================================

if __name__ == '__main__':
    if '--video' in sys.argv:
        generate_animation()
    else:
        generate_all_figures()
