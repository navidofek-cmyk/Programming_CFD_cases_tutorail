# Programming CFD Cases Tutorial

Educational 2D CFD mini-cases in modern C++.

The goal of this repository is to build readable numerical examples step by step:

- simple physics
- compact code
- no external dependencies
- standard library only
- educational structure rather than production complexity

## Cases Overview

| Case | Description | Method | Status |
|---|---|---|---|
| `01_cavity_case` | Lid-driven cavity | Finite differences + projection | available |
| `02_laminar_channel_obstacle` | Steady channel flow with square obstacle | Cell-centered FVM + SIMPLE-like pseudo-time | available |
| `03_channel_obstacle_unsteady` | Unsteady channel flow with square obstacle | Cell-centered transient solver + PISO-like corrections | available |
| `04_laminar_channel_obstacle_rebuild` | Rebuilt steady obstacle solver | Staggered MAC-style layout + SIMPLE/SOR | available |
| `05_laminar_channel_simple` | Simple laminar channel / Poiseuille case | Staggered channel solver | available |
| `06_exact_riemann_solver_1d_euler` | Exact 1D Euler reference case | Exact Riemann solution | available |
| `07_laval_nozzle_quasi_1d` | Quasi-1D compressible Laval nozzle | Explicit MacCormack + area source term | available |

## 01_cavity_case

2D lid-driven cavity flow:

- incompressible Navier-Stokes
- finite differences on a structured Cartesian grid
- projection / pressure-correction style approach
- explicit momentum update
- CSV export
- VTK export for ParaView

```bash
cd 01_cavity_case
./run.sh
```

## 02_laminar_channel_obstacle

2D laminar channel flow with a square obstacle:

- incompressible laminar Navier-Stokes
- cell-centered finite volume method
- structured Cartesian mesh
- SIMPLE-like pseudo-time steady iterations
- obstacle represented by a masked solid region
- CSV export

This case is intentionally simple and educational. It is useful for studying:

- pressure-velocity coupling
- obstacle masking on a Cartesian grid
- wake formation behind a bluff body
- practical convergence behavior of a simple segregated solver

```bash
cd 02_laminar_channel_obstacle
make
./channel_obstacle
```

## 03_channel_obstacle_unsteady

Transient companion to the steady obstacle case:

- incompressible laminar Navier-Stokes
- structured Cartesian mesh with masked square obstacle
- explicit predictor step with pressure correction
- PISO-like multiple corrections per time step
- CSV export for final fields
- VTK snapshots for ParaView

This case is meant for studying actual wake development in time instead of forcing the obstacle wake into a steady-state interpretation.

```bash
cd 03_channel_obstacle_unsteady
make
./channel_obstacle_unsteady
```

## 04_laminar_channel_obstacle_rebuild

Rebuilt steady channel-with-obstacle case with a cleaner pressure-velocity arrangement:

- staggered / MAC-style variable placement
- SIMPLE-like steady coupling
- structured Cartesian mesh
- CSV output

```bash
cd 04_laminar_channel_obstacle_rebuild
make
./cfd_channel
```

## 05_laminar_channel_simple

Simple laminar channel-flow case intended as a cleaner reference problem:

- structured channel geometry
- staggered-style formulation
- useful for comparing against expected Poiseuille-like behavior

```bash
cd 05_laminar_channel_simple
make
./channel_flow
```

## 06_exact_riemann_solver_1d_euler

Reference exact solution for the 1D compressible Euler Riemann problem:

- exact star-region solve for `p*` and `u*`
- wave-pattern sampling at user-selected time
- useful as a benchmark for later approximate compressible solvers

```bash
cd 06_exact_riemann_solver_1d_euler
make
./exact_riemann
```

## 07_laval_nozzle_quasi_1d

Educational compressible nozzle-flow case:

- quasi-1D Euler equations
- variable area duct
- explicit MacCormack predictor-corrector
- CSV output for `rho`, `u`, `p`, `T`, `M`

```bash
cd 07_laval_nozzle_quasi_1d
make
./nozzle_solver
```

## Notes

- These solvers are educational mini-codes, not production CFD tools.
- Numerical choices are intentionally simplified and commented honestly.
- The focus is readability, experimentation, and learning core CFD ideas.
