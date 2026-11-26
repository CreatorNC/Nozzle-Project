# Rocket Nozzle Design and Analysis

A 1D analytical model and a finite-volume (quasi-1D Euler) solver for a converging-diverging rocket nozzle.
Built with Python, NumPy and Matplotlib.

## What is in this folder

| File | Purpose |
|---|---|
| `nozzle_lib.py` | Isentropic relations (A/A*, p/p0, T/T0), nozzle design point, conical nozzle geometry |
| `rocket_nozzle.py` | Analytical performance model: thrust, Isp, flow profiles, thrust vs altitude, Isp vs expansion ratio |
| `validate_isentropic.py` | Checks the relations against standard gas-dynamics table values |
| `nozzle_cfd1d.py` | Quasi-1D Euler solver (HLL flux, RK2) with a mesh refinement study |
| `nozzle_profile.csv` | Nozzle contour (x, r in mm) for SolidWorks or for meshing |
| `results.json` | Numbers from the latest run |
| `*.png` | Plots produced by the scripts |

Run (Python 3, `pip install numpy matplotlib`):

```
python validate_isentropic.py
python rocket_nozzle.py
python nozzle_cfd1d.py
```

## Problem

Chamber 5 MPa and 3000 K, exhaust gas gamma = 1.2 and R = 350 J/kg/K, throat radius 50 mm.
Design the nozzle to be fully expanded at sea level (101325 Pa) and study its performance.

## Design point (from the analytical model)

| Quantity | Value |
|---|---|
| Exit Mach number | 3.03 |
| Area ratio Ae/At | 6.98 |
| Exit radius | 132.1 mm |
| Exit velocity | about 2454 m/s |
| Mass flow | about 24.9 kg/s |
| Thrust | about 61 kN |
| Isp (sea level) | about 250 s |

Conical nozzle: 30 degree convergent half-angle, 15 degree divergent half-angle, inlet radius 75 mm (1.5 x throat radius).
Convergent length 43.3 mm, divergent length 306.3 mm.

## Step 1: Validation of the relations

`validate_isentropic.py` compares the code with standard isentropic-flow table values for gamma = 1.4 at Mach 1.5, 2.0 and 3.0.
The worst error is 0.014 %, which is the rounding of the table itself. The round trip Mach -> A/A* -> Mach is exact to machine precision.

## Step 2: Finite-volume solver and mesh refinement

The quasi-1D Euler equations with the area source term are solved with a first-order finite-volume method, HLL fluxes and RK2 time marching.
The inlet uses stagnation conditions, the outlet uses supersonic extrapolation, and the flow starts from cold gas at rest and develops to steady state.
The converged solution is compared with the exact isentropic solution.

| Cells | Mach L2 error | Exit Mach error | Mass-flow error | Thrust error |
|---|---|---|---|---|
| 50 | 5.55 % | -3.71 % | +1.66 % | +1.11 % |
| 100 | 3.14 % | -1.96 % | +1.21 % | +0.93 % |
| 200 | 1.45 % | -0.94 % | +0.40 % | +0.24 % |
| 400 | 0.80 % | -0.50 % | +0.33 % | +0.25 % |

The Mach error falls roughly in proportion to the cell size (observed order about 0.8 to 1.1), as expected for a first-order scheme.
See `cfd1d_vs_exact.png` and `mesh_convergence.png`.

## Limits of this model (state these honestly when you present it)

- It is quasi-1D: no boundary layer, no wall friction, no heat loss, no two-dimensional effects such as flow divergence at the exit.
- Constant gamma and an ideal gas; the throat is a sharp corner.
- First-order accurate, so it smooths gradients near the throat.
- It is not a 2D axisymmetric CFD simulation. If you later run a 2D case in OpenFOAM or Fluent, `nozzle_profile.csv` gives the contour, and you can compare its centerline Mach number with these results.

## Things you should be able to explain

1. How the area-Mach relation is derived, and why the flow chokes at the throat.
2. Why the thrust is `mdot * ve + (pe - pa) * Ae`, and what over- and under-expansion mean.
3. How the HLL flux and the area source term work in the solver, and why the mesh study shows first-order behaviour.
