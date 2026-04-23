# Task: 2D Compressible RANS Solver for RAE 2822

Build a 2D compressible Navier-Stokes solver with Spalart-Allmaras (SA)
turbulence model in C++17 for the RAE 2822 airfoil at transonic conditions.

Reference case: NASA NPARC Case 4 — viscous, M∞ = 0.729, AoA = 2.31°,
Re = 6.5 × 10⁶. Experimental data: Cook, McDonald, Firmin (AGARD AR-138, 1979).

**Relation to case 11:** Start from the euler2d codebase. Keep the C-grid
generator (Mesh.hpp/cpp), Roe flux (Flux.cpp), BC framework, and IO writers.
Replace the explicit SSP-RK3 with LU-SGS implicit integration, add viscous
fluxes, and add the SA transport equation as a 5th equation.

---

## 1. Physical setup

Non-dimensional form: ρ∞ = 1, p∞ = 1/γ, |u∞| = M∞, γ = 1.4.

| Quantity          | Value          |
|-------------------|----------------|
| Mach number M∞    | 0.729          |
| Angle of attack   | 2.31°          |
| Reynolds number   | 6.5 × 10⁶      |
| γ                 | 1.4            |
| Pr (laminar)      | 0.72           |
| Pr_t (turbulent)  | 0.9            |
| Chord             | 1.0            |

Dynamic viscosity: μ = M∞ / Re  (non-dim: μ∞ = 1/Re).
Sutherland law is optional; constant μ is acceptable for this case.

Governing equations — compressible RANS in conservative form:

```
∂U/∂t + ∇·F_inv = ∇·F_vis + S

U    = [ρ, ρu, ρv, ρE, ρν̃]ᵀ

F_inv = standard Euler inviscid flux
F_vis = viscous + Reynolds stress flux (via Boussinesq hypothesis)
S    = [0, 0, 0, 0, S_SA]ᵀ  (SA source term)
```

Total effective viscosity: μ_eff = μ + μ_t,  μ_t = ρ ν̃ f_v1.

---

## 2. Mesh requirements

**Critical:** the RANS mesh must resolve the viscous sublayer.

Target: y⁺ ≈ 1 at the first cell centre.

Estimate for RAE 2822 at Re = 6.5×10⁶:
- Skin-friction coefficient: C_f ≈ 0.003
- Wall shear stress: τ_w = C_f · 0.5 · ρ∞ · u∞² ≈ 0.003 × 0.5 × 0.729² ≈ 8×10⁻⁴
- Friction velocity: u_τ = √(τ_w/ρ∞) ≈ 0.028
- Viscous length: δ_v = μ/ρ/u_τ = (1/Re) / u_τ ≈ 6.7×10⁻⁶
- First cell height for y⁺=1: Δy₁ ≈ δ_v = **~7×10⁻⁶ chord**

Recommended grid: **512 × 128** (nci × ncj, real cells).

```ini
ni_wrap    = 513    # = 2*n_wake + n_airfoil - 1
nj_normal  = 129
n_airfoil  = 257    # 2x finer than Euler case
n_wake     = 128
r_far      = 25.0
tanh_beta  = 7.5    # stronger wall clustering (was 3.5 in Euler case)
```

Verify after generation: h₁_avg < 1×10⁻⁵ chord, AR_wall < 5000.

---

## 3. Numerical methods

### 3.1 Inviscid flux — Roe + shock sensor MUSCL (reuse from case 11)

Same Roe flux with Harten entropy fix and pressure-based shock sensor.
Reconstruct all 5 variables (ρ, u, v, p, ν̃) with van Albada limiter.

### 3.2 Viscous flux — central difference with Green-Gauss gradients

For each face, compute the viscous flux:

```
F_vis · n = [0,
             (τ_xx·nx + τ_xy·ny),
             (τ_xy·nx + τ_yy·ny),
             (τ_xx·u + τ_xy·v)·nx + (τ_xy·u + τ_yy·v)·ny + q_x·nx + q_y·ny,
             (μ + μ_t/σ) · (∂ν̃/∂x·nx + ∂ν̃/∂y·ny)]
```

where:
```
τ_xx = μ_eff · (4/3·∂u/∂x - 2/3·∂v/∂y)
τ_yy = μ_eff · (4/3·∂v/∂y - 2/3·∂u/∂x)
τ_xy = μ_eff · (∂u/∂y + ∂v/∂x)
q    = -μ_eff/Pr_eff · ∂T/∂x   (Pr_eff blends Pr and Pr_t)
```

T = γ·p/ρ  (non-dimensional temperature).

**Gradient reconstruction (Green-Gauss):** for each cell (i,j), compute
∇φ = (1/A) Σ_faces φ_face · n · len. Use the face-average of the two
adjacent cell-centre values as φ_face.

### 3.3 Spalart-Allmaras turbulence model

Transport equation for ν̃ (modified turbulent kinematic viscosity):

```
∂(ρν̃)/∂t + ∇·(ρũν̃) = P - D + (1/σ)∇·[(μ + ρν̃)∇ν̃] + (cb2/σ)ρ|∇ν̃|²
```

Constants (Allmaras et al. 1994):
```
cb1 = 0.1355,  cb2 = 0.622,  σ = 2/3
cv1 = 7.1
cw1 = cb1/κ² + (1+cb2)/σ,  cw2 = 0.3,  cw3 = 2.0,  κ = 0.41
ct3 = 1.2,  ct4 = 0.5
```

Production:  `P = cb1 · ρ · S̃ · ν̃`
where S̃ = S + ν̃/(κ²d²)·f_v2,  S = |Ω| (vorticity magnitude),
f_v2 = 1 - χ/(1+χ·f_v1),  χ = ν̃/ν,  f_v1 = χ³/(χ³+cv1³).

Destruction: `D = cw1 · ρ · f_w · (ν̃/d)²`
f_w = g·[(1+cw3⁶)/(g⁶+cw3⁶)]^(1/6),  g = r + cw2·(r⁶-r),  r = ν̃/(S̃κ²d²).

Wall distance d: precompute after mesh generation as minimum distance
from each cell centre to the nearest j=0 wall node.

Boundary conditions for ν̃:
- Wall: ν̃ = 0
- Farfield (inflow): ν̃_∞ = 3 × ν = 3/Re
- Farfield (outflow): extrapolate from interior
- Wake cut: mirror (same as other variables)

### 3.4 Implicit time integration — LU-SGS

Lower-Upper Symmetric Gauss-Seidel. Replaces the explicit SSP-RK3.
Allows CFL = 10–100, achieving convergence in ~5× fewer iterations.

For steady state, the implicit update solves at each pseudo-time step:

```
(D + L) D⁻¹ (D + U) ΔU = -R(U)
```

where D = diagonal (local spectral radius), L = lower triangular sweep,
U = upper triangular sweep.

**Algorithm (one pseudo-time step):**

1. Compute R(U) — full inviscid + viscous residual.
2. Forward sweep (i,j increasing):
   ```
   (Ω/Δt + Λ_i) ΔU*_ij = -R_ij - Σ_{(i',j')∈L} A⁺_{i'/j'} ΔU*_{i',j'}
   ```
   Approximate off-diagonal terms with spectral radius:
   `A⁺ ΔU ≈ λ_face/2 · ΔU` (diagonal approximation).

3. Backward sweep (i,j decreasing), using ΔU* from forward:
   ```
   ΔU_ij = ΔU*_ij - (Ω/Δt + Λ_i)⁻¹ Σ_{(i',j')∈U} A⁻_{i'/j'} ΔU_{i',j'}
   ```

4. Update: `U ← U + ω · ΔU` (ω = 0.8–1.0 relaxation).

Diagonal: `Λ_i = 0.5 · Σ_faces (|λ_inv| + |λ_vis|) · len`
where `|λ_vis| = 2·(μ_eff/ρ) · len / area` (viscous spectral radius).

### 3.5 Wall boundary conditions (no-slip)

At j=0 (wall), ghost cell j=-1:
```
ρ_g  = ρ_0              (density extrapolated)
u_g  = -u_0             (no-slip: zero velocity at face → mirror)
v_g  = -v_0
T_g  =  T_0             (adiabatic wall: ∂T/∂n = 0 → mirror)
p_g  =  p_0             (zero pressure gradient at wall)
ν̃_g  = -ν̃_0            (ν̃=0 at wall → mirror)
```

Adiabatic wall (no heat flux).  For isothermal wall, set T_wall = T_stag.

### 3.6 Farfield BC

Same Riemann characteristic BC as Euler case, extended to 5 variables.
For ν̃: inflow ghost = 3/Re, outflow ghost = extrapolated.

---

## 4. Project structure

```
rans2d/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── Airfoil.hpp     (copy from euler2d)
│   ├── Mesh.hpp        (copy + add tanh_beta param)
│   ├── Solver.hpp      (new: 5-variable, LU-SGS)
│   ├── Flux.hpp        (copy Roe flux; add viscous flux)
│   ├── BC.hpp
│   ├── IO.hpp
│   └── SpalartAllmaras.hpp
├── src/
│   ├── main.cpp
│   ├── Mesh.cpp        (copy + update)
│   ├── Solver.cpp
│   ├── Flux.cpp        (copy + add viscous)
│   ├── BC.cpp
│   ├── IO.cpp
│   └── SpalartAllmaras.cpp
└── config/
    └── rae2822.ini
```

---

## 5. Config file

```ini
# Flow
mach        = 0.729
aoa_deg     = 2.31
gamma       = 1.4
reynolds    = 6.5e6
prandtl     = 0.72
prandtl_t   = 0.9

# Grid
ni_wrap     = 513
nj_normal   = 129
n_airfoil   = 257
n_wake      = 128
r_far       = 25.0
tanh_beta   = 7.5

# Solver
cfl         = 20.0          # LU-SGS allows large CFL
max_iter    = 5000
residual_drop = 1.0e-6
scheme_order = 2
warmup_iters = 500
output_interval = 200
```

---

## 6. Acceptance criteria

1. Build without warnings (`-Wall -Wextra`).
2. `./rans2d config/rae2822.ini` produces `output/field.vts`, `output/cp.dat`,
   `output/convergence.dat`.
3. Density residual drops ≥ 5 orders.
4. Cp distribution matches AGARD AR-138 experiment within ±0.05 Cp
   (viscous shock position closer to x/c ≈ 0.55 vs. Euler's 0.61).
5. Skin-friction coefficient C_f > 0 everywhere (no separation predicted
   for this mild condition with SA model).
6. No NaNs in the field; ν̃ ≥ 0 everywhere.

---

## 7. Execution plan

### Stage 1 — Mesh
- Generate fine RANS C-grid (512×128 cells, h₁ ≈ 7×10⁻⁶).
- Verify y⁺ ≈ 1 estimate.
- Checkpoint: `mesh.vts` opens in ParaView, wall spacing is correct.

### Stage 2 — LU-SGS on Euler equations
- Implement LU-SGS implicit solver for the 4-variable Euler system.
- Verify: same Cp as euler2d case 11, but converges in ~500–1000 iters.
- Checkpoint: Cp matches, residual 5 orders in < 2000 iters.

### Stage 3 — Viscous fluxes (laminar NS)
- Add Green-Gauss gradient reconstruction.
- Add viscous flux contributions to the residual.
- No turbulence model yet (laminar).
- Checkpoint: C_f distribution plausible, boundary layer thickness visible
  in wall-normal profiles.

### Stage 4 — Spalart-Allmaras
- Add 5th equation (ρν̃) to the system.
- Implement SA production, destruction, diffusion, wall distance.
- Checkpoint: ν̃ > 0 everywhere in the boundary layer,
  Cp shifts toward experimental data, C_f ≥ 0.

### Stage 5 — Output and comparison
- Plot Cp vs AGARD experiment.
- Plot C_f distribution.
- Extract drag (pressure + friction) and lift coefficients.
- Compare CL with viscous experiment (~0.743).

---

## 8. Rules

- C++17, no external libraries.
- Reuse Airfoil.hpp, Mesh.cpp, Roe flux, IO writers from euler2d.
- Commit after each stage.
- Ask before deviating from this spec.
