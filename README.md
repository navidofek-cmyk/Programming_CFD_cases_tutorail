# Programming CFD Cases Tutorial

Small educational CFD cases in modern C++.

The goal of this repository is to build readable numerical examples step by step:

- simple physics
- compact code
- no external dependencies
- standard library only
- educational structure rather than production complexity

## Cases

### 01_cavity_case

2D lid-driven cavity flow:

- incompressible Navier-Stokes
- finite differences on a structured Cartesian grid
- projection / pressure-correction style approach
- explicit momentum update
- CSV export
- VTK export for ParaView

Good as a first CFD teaching example because the geometry and boundary conditions are simple.

### 02_laminar_chanell_obstacle

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

## Build and Run

Each case has its own local instructions inside its folder.

Typical usage:

```bash
cd 01_cavity_case
./run.sh
```

or:

```bash
cd 02_laminar_chanell_obstacle
make
./channel_obstacle
```

## Notes

- These solvers are educational mini-codes, not production CFD tools.
- Numerical choices are intentionally simplified and commented honestly.
- The focus is readability, experimentation, and learning core CFD ideas.
