# 2D Euler solver — RAE 2822 transonic airfoil

C++17 cell-centred finite-volume solver for the 2D compressible Euler equations on a structured C-grid.  
Reference case: M∞ = 0.729, AoA = 2.31°, γ = 1.4 (inviscid, non-dimensional).

---

## Numerical methods

| Component | Method |
|-----------|--------|
| Grid | C-topology structured, TFI + tanh wall clustering |
| Spatial discretisation | Cell-centred FV, 2nd order MUSCL + van Albada limiter |
| Numerical flux | Roe scheme with Harten entropy fix (δ = 0.1 â) |
| Time integration | SSP-RK3, local time stepping |
| Boundary conditions | Inviscid wall (velocity reflection), Riemann characteristic farfield, wake-cut mirror |
| Convergence strategy | First-order warmup → shock-sensor blended MUSCL |

---

## Result

- Transonic shock on upper surface at **x/c ≈ 0.61**
- Cp_min ≈ −1.31 on the suction side
- Residual drops **~5 orders of magnitude** from freestream initialisation

---

## MUSCL limit-cycle problem and fix

### What went wrong

When the solver switched from first-order (warmup) to full MUSCL reconstruction,
the residual quickly plateaued at ~3 × 10⁻³ and refused to drop further — a
**limit cycle** rather than convergence.

**Root cause:** the transonic normal shock spans 1–2 cells.  
The van Albada limiter is a smooth, symmetric function of the local gradient ratio.
At the shock, the limiter allows partial reconstruction (≈ 25–50 % of the slope).
Every RK iteration:

1. MUSCL reconstructs slightly different left/right states at the shock face.
2. The Roe flux drives the shock one half-cell downstream.
3. The next iteration drives it one half-cell back.
4. The oscillation is self-sustaining because the explicit time stepping has no
   mechanism to damp it — the fixed-point of the MUSCL scheme *at the shock*
   is not a steady-state but a 2-cycle.

This is a known failure mode of explicit MUSCL on transonic/supersonic shocks.
It has nothing to do with stability (the scheme is stable) and cannot be cured by
reducing the CFL.

### Fix: pressure-based shock sensor

For every face, compute a normalised pressure jump across the 4-cell MUSCL stencil:

```
s = max( |pA−p0|/(pA+p0),  |p0−p1|/(p0+p1),  |p1−pB|/(p1+pB) )
```

Then scale the reconstruction factor locally:

```
alpha(s) = max(0,  1 − 3·s)
```

| Location | Δp/p | alpha | Effective order |
|----------|------|-------|-----------------|
| Transonic shock | ≈ 0.38 | 0 | 1st (no oscillation) |
| Smooth supersonic region | ≈ 0.005 | 0.98 | ≈ 2nd (full accuracy) |
| Boundary layer / LE | ≈ 0.05 | 0.85 | ≈ 2nd |

At the shock the reconstruction is suppressed to first order, eliminating the
oscillation source.  Everywhere else the scheme remains second order.  
The Cp distribution is unchanged because the shock is already sharp at first order
(Rankine–Hugoniot is satisfied cell-by-cell); only the smooth upstream/downstream
regions benefit from the second-order reconstruction.

**Outcome:** residual drops from 3 × 10⁻³ plateau to ~3 × 10⁻⁵ — five full orders
from initial freestream, vs. three orders without the sensor.

---

## Usage

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/euler2d config/rae2822.ini
python3 plot_cp.py
python3 plot_field.py
python3 plot_convergence.py
```

Output files in `output/`:

| File | Content |
|------|---------|
| `field.vts` | VTK StructuredGrid — full flow field (ParaView) |
| `mesh.vts` | VTK StructuredGrid — mesh quality (area, AR) |
| `surface.vtp` | VTK PolyData — airfoil surface with Cp |
| `cp.dat` | ASCII (x/c, Cp) for quick plotting |
| `convergence.dat` | Iteration history |

---

## Key config parameters (`config/rae2822.ini`)

```ini
warmup_iters     = 5000   # first-order iterations before MUSCL
cfl              = 1.5    # CFL during warmup
cfl_muscl        = 0.8    # CFL after MUSCL switch
muscl_ramp_iters = 0      # 0 = instant switch; >0 = gradual ramp
residual_drop    = 1e-5   # stop when res/res0 < this
```
