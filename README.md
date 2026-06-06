# Qian Xuesen Extended Trajectory Optimization

## Project Overview

Extend Qian Xuesen's patched conic trajectory design (*Introduction to Interstellar Navigation*, Ch.5-6) using modern numerical methods and JPL Horizons ephemeris data.

**Mission:** Earth → Moon gravity assist → heliocentric transfer to near-Sun perihelion → return to Earth. Find the optimal 2026 launch date minimizing total Δv.

## Quick Start

### Environment
```bash
conda create -n final_project python=3.10
conda activate final_project
pip install numpy scipy matplotlib astroquery astropy
```

### Run All
```bash
make all        # Generate figures + report PDF
```

### Run Individual Steps
```bash
# Step 1: Verify analytical solution
python src/patched_conic.py

# Step 2: Run full-year optimization (coarse scan, ~15s)
python src/optimize_launch.py

# Step 3: Run with fine scan + refinement (~3min)
python src/optimize_launch.py --fine

# Step 4: Sensitivity analysis
python src/sensitivity.py

# Step 5: Generate all figures
python src/visualize.py

# Step 6: Compile report
make pdf
```

## Project Structure

```
qian-xuesen-trajectory-optimization/
├── .gitignore
├── Makefile                  # make all / make pdf / make clean
├── README.md                 # This file
├── AI-Agent.md               # AI tool usage log
├── report.tex                # Final report (XeLaTeX)
├── report.pdf                # Compiled report
├── src/
│   ├── physical_constants.py # Unified constants & constraints
│   ├── patched_conic.py      # M1: Patched conic analytical solution
│   ├── nbody.py              # M2: Velocity-Verlet 4-body integrator
│   ├── horizons_validate.py  # M3: JPL Horizons validation
│   ├── lunar_flyby.py        # M4: Vector-based lunar flyby model
│   ├── single_trajectory.py  # M5: Full trajectory solver
│   ├── optimize_launch.py    # M6: 3-phase launch window optimizer
│   ├── sensitivity.py        # M7: Parameter sensitivity analysis
│   └── visualize.py          # M8: 5 types of visualizations
├── data/
│   ├── scan_2026_coarse.json # Coarse scan: 365 days × 8 rp values
│   ├── scan_2026_fine.json   # Fine scan: top-10 days × (rm, rp) grid
│   ├── scan_2026_optimal.json# Refined optimal solution
│   ├── dv_vs_day.json        # Δv vs launch day curve data
│   └── sensitivity.json      # Sensitivity analysis results
└── figures/
    ├── trajectory.png        # Optimal trajectory 2D plot
    ├── trajectory.mp4        # 45s trajectory animation (O5)
    ├── energy_conservation.png
    ├── dv_vs_day.png         # Δv vs launch day
    ├── horizons_errors.png   # N-body vs ephemeris errors
    └── sensitivity_heatmap.png
```

## Module Status

| Module | Description | Status | Author |
|---|---|---|---|
| M1 | Patched Conic Analytical | ✅ | DeepSeek |
| M2 | N-body Integrator | ✅ | DeepSeek |
| M3 | Horizons Validation | ✅ | DeepSeek |
| M4 | Lunar Gravity Assist | ✅ | Claude Opus |
| M5 | Single Trajectory Solver | ✅ | Claude Opus |
| M6 | Launch Window Optimization | ✅ | DeepSeek (framework) + Claude Opus (analysis) |
| M7 | Sensitivity Analysis | ✅ | Claude Opus |
| M8 | Visualization | ✅ | DeepSeek |
| Report | LaTeX | ✅ | Claude Opus |

## Key Results

### Optimal Trajectory (2026)
- **Launch day:** Day 107 (April 18, 2026)
- **Moon flyby distance:** 2245 km (trailing side)
- **Perihelion:** 0.400 AU
- **Δv_total:** 26.177 km/s
- **Savings vs direct:** 0.029 km/s (0.11%)

### M1 Validation (rp = 0.2 AU, classic Qian case)
- Heliocentric transfer: a = 0.600 AU, e = 0.667
- Δv_total (direct) = 33.21 km/s, matches Qian §6.3

### Sensitivity (at optimal point)
- ∂Δv/∂rp ≈ −30 km/s/AU (dominant)
- ∂Δv/∂rm ≈ −0.002 km/s/1000km (negligible)
- ∂Δv/∂t₀ ≈ 0.032 km/s/day
- Moon phase error > ±2° negates the 0.03 km/s savings

### Physical Conclusion
Lunar gravity assist provides < 0.5% Δv savings for this mission profile. The Moon's mass is too small relative to the Sun, and the high-energy trajectory (v∞ ≫ v_escape,Moon) limits the flyby deflection angle to < 3°.

## References

- Qian Xuesen, *Introduction to Interstellar Navigation*, Ch.5-6 (2008 ed.)
- Vallado D.A., *Fundamentals of Astrodynamics and Applications*, 4th ed., §12.4
- Curtis H.D., *Orbital Mechanics for Engineering Students*, 3rd ed., §8.10
- JPL Horizons: https://ssd.jpl.nasa.gov/horizons/manual.html
