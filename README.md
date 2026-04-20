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
| `02_laminar_chanell_obstacle` | Steady channel flow with square obstacle | Cell-centered FVM + SIMPLE-like pseudo-time | available |
| `03_channel_obstacle_unsteady` | Unsteady channel flow with square obstacle | Cell-centered transient solver + PISO-like corrections | available |

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

## 02_laminar_chanell_obstacle

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
cd 02_laminar_chanell_obstacle
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

## Notes

- These solvers are educational mini-codes, not production CFD tools.
- Numerical choices are intentionally simplified and commented honestly.
- The focus is readability, experimentation, and learning core CFD ideas.
