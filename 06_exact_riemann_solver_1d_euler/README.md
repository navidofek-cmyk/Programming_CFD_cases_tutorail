# 06_exact_riemann_solver_1d_euler

Educational exact Riemann solver for the 1D compressible Euler equations.

## Purpose

This case solves a classic 1D shock tube problem by evaluating the exact
self-similar Riemann solution, not by time marching a finite-volume scheme.

That means:

- left and right initial states are separated by a diaphragm
- the solver computes the star-region pressure `p*` and velocity `u*`
- the final solution at time `t` is sampled analytically from the wave pattern

## Default Problem

The default setup is the Sod shock tube:

- left:  `rho = 1.0`, `u = 0.0`, `p = 1.0`
- right: `rho = 0.125`, `u = 0.0`, `p = 0.1`
- `gamma = 1.4`

## Numerical Idea

The implementation follows the standard exact Riemann procedure:

1. solve the nonlinear equation for `p*` using Newton iteration
2. recover `u*`
3. identify whether each side produces a shock or a rarefaction
4. sample the self-similar solution as a function of `s = (x - x0) / t`

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/06_exact_riemann_solver_1d_euler
make
```

## Run

```bash
./exact_riemann
```

or:

```bash
make run
```

## Output

The solver writes:

- `solution.csv`

with columns:

- `x`
- `rho`
- `u`
- `p`
- `e`

where `e` is the specific internal energy.

## Notes

- This is an exact Riemann solution sampler, not a full CFD solver.
- It is ideal as a reference solution for later approximate solvers:
  - Rusanov
  - Roe
  - HLL / HLLC
- It is also useful for checking wave structure:
  - rarefaction
  - contact discontinuity
  - shock
