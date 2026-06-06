"""
nbody.py  M2: N 

Velocity-Verlet ,  Sun-Earth-Moon-Rocket 

:
- Velocity-Verlet  (2),  3600s,  60s
-  (Rocket, body)
- Energy conservation ( 1000 )
- 2-body benchmarkVerification: =4 , 1 yrrel_error  10

:  km,  s,  km/s
: J2000  (2D: x-y )

:
    from nbody import NBodySimulation, run_two_body_benchmark
    sim = NBodySimulation()
    sim.set_ics_from_horizons(t0_jd)  #  Horizons IC
    history = sim.integrate(t_span=365*86400, h=3600.0)
"""

import numpy as np
from physical_constants import (
    MU_SUN, MU_EARTH, MU_MOON,
    AU, R_EARTH, R_MOON, R_EARTH_ORBIT, R_MOON_ORBIT,
    SECONDS_PER_DAY, SECONDS_PER_YEAR,
    MU_NORM, C_ENERGY_TOL,
)


# ============================================================
#  1. 
# ============================================================

def compute_accelerations(positions, masses, mu_values):
    """
    body ()

    Rocket (, mass  0): body, 

    Parameters
    ----------
    positions : ndarray (N, 2)
        body [x, y] (km)
    masses : ndarray (N,)
        Rocket 0 ()
    mu_values : ndarray (N,)
          = G*m [km/s]Rocket  = 0

    Returns
    -------
    ndarray (N, 2): body [ax, ay] (km/s)
    """
    N = len(positions)
    acc = np.zeros_like(positions)

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            r_vec = positions[j] - positions[i]
            r = np.linalg.norm(r_vec)

            # 
            if r < 1e-6:
                continue

            # _j / r * 
            acc[i] += mu_values[j] * r_vec / (r**3)

    return acc


def compute_accelerations_vectorized(positions, mu_values):
    """
     (body)

    Parameters
    ----------
    positions : ndarray (N, 2)
    mu_values : ndarray (N,)

    Returns
    -------
    ndarray (N, 2)
    """
    N = len(positions)
    acc = np.zeros_like(positions)

    for i in range(N):
        r_vec = positions - positions[i]  # (N, 2)
        r = np.linalg.norm(r_vec, axis=1)  # (N,)
        r = np.where(r < 1e-6, np.inf, r)  # 
        acc[i] = np.sum(
            mu_values[:, np.newaxis] * r_vec / (r[:, np.newaxis]**3),
            axis=0
        )

    return acc


# ============================================================
#  2. Velocity-Verlet 
# ============================================================

class VelocityVerlet:
    """
    Velocity-Verlet  (, 2 , )

    Algorithm:
        x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt
        a(t+dt) = compute_accel(x(t+dt))
        v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    """

    def __init__(self, compute_accel_fn):
        """
        Parameters
        ----------
        compute_accel_fn : callable(pos, mu_values) -> acc
            
        """
        self.compute_accel = compute_accel_fn

    def step(self, pos, vel, mu_values, dt, acc_current=None):
        """
        

        Parameters
        ----------
        pos, vel : ndarray (N, 2)
            
        mu_values : ndarray (N,)
            
        dt : float
             [s]
        acc_current : ndarray or None
             ()

        Returns
        -------
        pos_new, vel_new, acc_new
        """
        if acc_current is None:
            acc_current = self.compute_accel(pos, mu_values)

        # 
        pos_new = pos + vel * dt + 0.5 * acc_current * dt**2

        # 
        acc_new = self.compute_accel(pos_new, mu_values)

        #  ()
        vel_new = vel + 0.5 * (acc_current + acc_new) * dt

        return pos_new, vel_new, acc_new


# ============================================================
#  3. N 
# ============================================================

class NBodySimulation:
    """
    N 

    body: [Sun, Earth, Moon, Rocket]
    Rocket : mass=0, =0.
    """

    # body
    BODY_NAMES = ['Sun', 'Earth', 'Moon', 'Rocket']

    def __init__(self):
        self.integrator = VelocityVerlet(
            lambda pos, mu: compute_accelerations(pos, None, mu)
        )
        self.positions = None   # (4, 2) [km]
        self.velocities = None  # (4, 2) [km/s]
        self.mu_values = None   # (4,) [km/s]
        self.t_current = 0.0    # [s]
        self.history = []       # Energy

    def set_bodies(self, pos_sun, vel_sun,
                    pos_earth, vel_earth,
                    pos_moon, vel_moon,
                    pos_rocket=None, vel_rocket=None):
        """
        body

         (2,)  (3,) ,  km, km/s
        Rocket, Earth ()
        """
        N = 4
        self.positions = np.zeros((N, 2))
        self.velocities = np.zeros((N, 2))
        self.mu_values = np.zeros(N)

        # Sun
        self.positions[0] = pos_sun[:2]
        self.velocities[0] = vel_sun[:2]
        self.mu_values[0] = MU_SUN

        # Earth
        self.positions[1] = pos_earth[:2]
        self.velocities[1] = vel_earth[:2]
        self.mu_values[1] = MU_EARTH

        # Moon
        self.positions[2] = pos_moon[:2]
        self.velocities[2] = vel_moon[:2]
        self.mu_values[2] = MU_MOON

        # Rocket ()
        if pos_rocket is None:
            self.positions[3] = pos_earth[:2].copy()  # Earth
        else:
            self.positions[3] = pos_rocket[:2]
        if vel_rocket is None:
            self.velocities[3] = vel_earth[:2].copy()
        else:
            self.velocities[3] = vel_rocket[:2]
        self.mu_values[3] = 0.0  # 

        self.t_current = 0.0

    def set_ics(self, pos, vel, mu):
        """IC

        pos: (N,2), vel: (N,2), mu: (N,)
        """
        N = len(pos)
        self.positions = np.array(pos, dtype=float)
        self.velocities = np.array(vel, dtype=float)
        self.mu_values = np.array(mu, dtype=float)
        self.t_current = 0.0

    def compute_total_energy(self):
        """
        Energy ( + )

        E =  0.5*m_i*v_i - _{i<j} G*m_i*m_j / r_ij

        Rocket (m=0, =0), 
        body + RocketEnergybodyOK /G 

        ,    Gm:
        E/m_i = 0.5 *  v_i - _{i<j} _j / r_ij  (Energy)
        """
        N = len(self.positions)
        E_kin = 0.0
        E_pot = 0.0

        for i in range(N):
            # : 0.5 * v_i
            # :  (mu=0), 
            E_kin += 0.5 * np.sum(self.velocities[i]**2)

            for j in range(i + 1, N):
                r = np.linalg.norm(self.positions[i] - self.positions[j])
                if r < 1e-6:
                    continue
                # : -_i * _j / r ...
                # body i  j: -_i * (1) / r
                # Simplified - _i / r  ()
                if self.mu_values[j] > 0:
                    E_pot -= self.mu_values[j] / r
                if self.mu_values[i] > 0:
                    E_pot -= self.mu_values[i] / r

        return E_kin + E_pot

    def compute_rocket_energy(self):
        """
        RocketEnergy (, )
        Energy conservationcheck (C5)
        """
        v_sq = np.sum(self.velocities[3]**2)  # v
        E_kin = 0.5 * v_sq

        E_pot = 0.0
        for i in range(3):  # Sun, Earth, Moon
            r = np.linalg.norm(self.positions[3] - self.positions[i])
            if r > 1e-6:
                E_pot -= self.mu_values[i] / r

        return E_kin + E_pot

    def integrate(self, t_span, h=3600.0, h_fine=60.0,
                   fine_body_idx=None, fine_radius=50000.0,
                   energy_check_interval=1000,
                   callback=None):
        """
         N 

        Parameters
        ----------
        t_span : float
            Time [s]
        h : float
             [s] ( 3600s)
        h_fine : float
             [s] ( 60s, body)
        fine_body_idx : int or None
            body (: Moon=2)
        fine_radius : float
            radius [km] (body <  h_fine)
        energy_check_interval : int
            checkEnergy
        callback : callable(t, pos, vel) or None
             ()

        Returns
        -------
        dict: {
            't_history': [...],
            'energy_history': [...],
            'pos_history': [...]  ( callback )
        }
        """
        if fine_body_idx is None:
            fine_body_idx = 2  # : Moon

        N = len(self.positions)
        steps_total = int(np.ceil(t_span / h))

        # 
        acc = self.integrator.compute_accel(self.positions, self.mu_values)

        # Energy
        energy_history = [(0.0, self.compute_rocket_energy())]
        t_history = [0.0]

        t = 0.0
        step_count = 0

        while t < t_span:
            # : Rocketbody
            if fine_body_idx is not None:
                r_rocket_to_body = np.linalg.norm(
                    self.positions[3] - self.positions[fine_body_idx]
                )
                dt = h_fine if r_rocket_to_body < fine_radius else h
            else:
                dt = h

            # 
            if t + dt > t_span:
                dt = t_span - t

            # 
            self.positions, self.velocities, acc = self.integrator.step(
                self.positions, self.velocities, self.mu_values, dt, acc
            )

            t += dt
            step_count += 1
            self.t_current = t

            # Energy
            if step_count % energy_check_interval == 0:
                E = self.compute_rocket_energy()
                energy_history.append((t, E))
                t_history.append(t)

            #  ()
            if callback is not None:
                callback(t, self.positions.copy(), self.velocities.copy())

        return {
            't_history': np.array(t_history),
            'energy_history': np.array(energy_history),
            'n_steps': step_count,
        }

    def apply_impulse(self, body_idx, delta_v):
        """
        body ()

        Parameters
        ----------
        body_idx : int
            body (3 = Rocket)
        delta_v : ndarray (2,)
             [km/s]
        """
        self.velocities[body_idx] += np.array(delta_v[:2])

    def get_rocket_state(self):
        """Rocket"""
        return self.positions[3].copy(), self.velocities[3].copy()


# ============================================================
#  4. 2-body benchmarkVerification
# ============================================================

def run_two_body_benchmark():
    """
    Verification Velocity-Verlet 

    :  = 4,  m=1,  = AU,  = year/(2)  58.13 .
    :  r=1, v=2, T=1 ( 2 TU = 1 year).

    Verification:  1 yr (2 TU), rel_error  10.
    """
    print("=" * 60)
    print("  M2 Verification: Velocity-Verlet 2-body benchmark")
    print("=" * 60)

    # 
    mu = MU_NORM  # 4

    # IC:  r=1, v=2
    r0 = 1.0
    v0 = 2.0 * np.pi

    pos = np.array([
        [0.0, 0.0],    # body (Sun, )
        [r0, 0.0],     # ,  (1, 0) 
    ], dtype=float)

    vel = np.array([
        [0.0, 0.0],
        [0.0, v0],     #  +y 
    ], dtype=float)

    mu_vals = np.array([mu, 0.0])  # 

    sim = NBodySimulation()
    sim.set_ics(pos, vel, mu_vals)

    h_norm = 0.001  # fine step for 1e-4 accuracy
    T_norm = 2.0 * np.pi  # 1 yr = 2 TU

    acc_fn = lambda p, m: compute_accelerations_vectorized(p, m)

    integrator = VelocityVerlet(acc_fn)
    acc = integrator.compute_accel(pos, mu_vals)

    t = 0.0
    positions = [(t, pos.copy())]
    pos_cur = pos.copy()
    vel_cur = vel.copy()

    while t < T_norm:
        dt = h_norm
        if t + dt > T_norm:
            dt = T_norm - t

        pos_cur, vel_cur, acc = integrator.step(pos_cur, vel_cur, mu_vals, dt, acc)
        t += dt

        if t % 0.5 < h_norm:  #  0.5 TU 
            positions.append((t, pos_cur.copy()))

    # 
    r_final = np.linalg.norm(pos_cur[1])
    r_expected = r0  #  r=1

    # check: , r 
    r_max = max(np.linalg.norm(p[1]) for _, p in positions)
    r_min = min(np.linalg.norm(p[1]) for _, p in positions)

    error_rms = np.std([np.linalg.norm(p[1]) - r0 for _, p in positions])
    error_rel = (r_max - r_min) / r0

    print(f"  Step h = {h_norm:.3f} TU")
    print(f"  Time T = {T_norm:.3f} TU (1 yr)")
    print(f"  radius r = {r0:.6f}")
    print(f"  radius r = {r_final:.6f}")
    print(f"  radiusrange: [{r_min:.8f}, {r_max:.8f}]")
    print(f"  rel_error = {error_rel:.2e}")
    print(f"  RMS  = {error_rms:.2e}")

    passed = error_rel <= 1e-4
    print(f"\n  Result: {'PASS' if passed else 'FAIL'} (threshold 1e-4)")

    return {
        'error_rel': error_rel,
        'error_rms': error_rms,
        'passed': passed,
    }


# ============================================================
#  5. Sun-Earth-Moon 3-bodyVerification
# ============================================================

def run_three_body_verification(t_span_days=365, h=3600.0):
    """
    IC Sun-Earth-Moon , Verification

    IC: Simplified ( JPL )
     Horizons Verification M3 
    """
    print("\n" + "=" * 60)
    print("  M2: Sun-Earth-Moon 3-body (Simplified IC)")
    print("=" * 60)

    # SimplifiedIC
    # Sun
    pos_sun = np.array([0.0, 0.0])
    vel_sun = np.array([0.0, 0.0])

    # Earth (AU, 0),  +y
    pos_earth = np.array([AU, 0.0])
    vel_earth = np.array([0.0, np.sqrt(MU_SUN / AU)])

    # MoonEarth + (R_moon_orbit, 0),  = Earth + 
    pos_moon = np.array([AU + R_EARTH_ORBIT, 0.0]) if False else np.array([AU, R_MOON_ORBIT])
    # Simplified: MoonEarth x 
    pos_moon = np.array([AU + 384400.0, 0.0])  # 
    v_moon_earth = np.sqrt(MU_EARTH / 384400.0)
    vel_moon = np.array([0.0, np.sqrt(MU_SUN / AU) + v_moon_earth])

    # Rocket: Earth
    pos_rocket = pos_earth.copy()
    vel_rocket = vel_earth.copy()

    sim = NBodySimulation()
    sim.set_bodies(pos_sun, vel_sun, pos_earth, vel_earth, pos_moon, vel_moon,
                    pos_rocket, vel_rocket)

    t_span = t_span_days * SECONDS_PER_DAY

    # 
    trajectory = []

    def record_traj(t, pos, vel):
        trajectory.append((t, pos.copy(), vel.copy()))

    result = sim.integrate(
        t_span, h=h, h_fine=60.0,
        fine_body_idx=2,  # Moon
        energy_check_interval=100,
        callback=record_traj,
    )

    # Energy conservationcheck
    E_init = result['energy_history'][0, 1]
    E_final = result['energy_history'][-1, 1]
    E_rel_error = abs((E_final - E_init) / E_init) if abs(E_init) > 1e-12 else 0.0

    print(f"  Time: {t_span_days} ")
    print(f"  Steps: {result['n_steps']}")
    print(f"  Initial rocketEnergy: {E_init:.6f} km/s")
    print(f"  Final rocketEnergy: {E_final:.6f} km/s")
    print(f"  Energyrel_error: {E_rel_error:.2e}")
    print(f"  Energycheckpoints: {len(result['energy_history'])}")

    energy_ok = E_rel_error <= C_ENERGY_TOL
    print(f"  C5 Energy conservation: {'PASS' if energy_ok else 'FAIL'} (threshold {C_ENERGY_TOL})")

    return {
        'energy_rel_error': E_rel_error,
        'energy_ok': energy_ok,
        'n_steps': result['n_steps'],
        'trajectory': trajectory,
    }


# ============================================================
#  Main
# ============================================================

if __name__ == '__main__':
    # 2-body benchmark
    bench = run_two_body_benchmark()

    # 3-bodySimplifiedVerification
    three_body = run_three_body_verification(t_span_days=30, h=3600.0)

    print()
    print("M2 .")
