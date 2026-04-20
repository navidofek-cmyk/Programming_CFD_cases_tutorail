# 03_channel_obstacle_unsteady

Small educational transient CFD case for 2D incompressible flow in a channel with a square obstacle.

## Purpose

This case is the unsteady companion to the steady obstacle case:

- `02_laminar_chanell_obstacle` stays as a steady SIMPLE-like solver
- `03_channel_obstacle_unsteady` treats the wake as a real time-dependent problem

The code uses:

- C++20
- standard library only
- structured Cartesian mesh
- cell-centered fields
- obstacle masking on the grid
- PISO-like pressure-velocity coupling

## Geometry

- channel length `Lx = 32H`
- channel height `Ly = 8H`
- square obstacle from `x = 8H` to `9H`
- square obstacle from `y = 3.5H` to `4.5H`
- `H = 1`

## Physics

- incompressible laminar Navier-Stokes
- constant density
- Reynolds number `Re = 60`
- inlet velocity `U_in = 1`
- viscosity computed from `nu = U_in * H / Re`

## Numerics

This is an educational PISO-like solver, not a full production PISO implementation.

Per time step:

1. build a provisional velocity from explicit convection and diffusion plus the old pressure gradient
2. solve a pressure Poisson equation from the divergence of the provisional velocity
3. correct pressure and velocity
4. repeat the pressure correction sweep a small number of times

Main simplifications:

- collocated cell-centered layout
- explicit momentum update
- simple pressure Poisson solve with fixed iteration count
- obstacle represented by a stair-step solid mask
- no higher-order convection scheme

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/03_channel_obstacle_unsteady
make
```

## Run

```bash
./channel_obstacle_unsteady
```

or:

```bash
make run
```

## Output

Final fields:

- `u.csv`
- `v.csv`
- `p.csv`

VTK snapshots for ParaView:

- `output/flow_00100.vtk`
- `output/flow_00200.vtk`
- ...

Console output includes:

- time step
- physical time
- change in `u`
- change in `v`
- pressure-correction residual
- divergence norm
- maximum velocity magnitude

## ParaView

Open any snapshot, for example:

```bash
paraview output/flow_00600.vtk
```

Useful views:

- pressure contours
- velocity glyphs
- stream tracer

## Notes

This case is intended for learning:

- wake development in time
- pressure-velocity coupling for incompressible flow
- how obstacle masking works on a Cartesian grid

It is deliberately simple and honest about the numerical approximations.
