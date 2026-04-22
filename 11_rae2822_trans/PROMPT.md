# Task: 2D Euler (inviscid) Finite-Volume Solver for RAE 2822

You are going to build a working 2D compressible Euler solver in C++17
that computes steady transonic flow around the RAE 2822 airfoil. The
reference validation case is NASA Glenn's NPARC archive:
https://www.grc.nasa.gov/www/wind/valid/raetaf/raetaf04/raetaf04.html

**Important scope clarification:** The NASA page describes a *viscous*
(RANS) calculation. We are doing the **inviscid (Euler) version only**.
Ignore Reynolds number, turbulence models, viscous walls, and y+
discussion from that page. We keep Mach, AoA, and the airfoil geometry.

---

## 1. Physical setup

Freestream (non-dimensional form — use `ρ∞ = 1`, `p∞ = 1/γ`,
`|u∞| = M∞`, γ = 1.4):

| Quantity        | Value   |
|-----------------|---------|
| Mach number M∞  | 0.729   |
| Angle of attack | 2.31°   |
| γ               | 1.4     |
| Chord           | 1.0     |

Governing equations: 2D compressible Euler in conservative form

```
∂U/∂t + ∂F(U)/∂x + ∂G(U)/∂y = 0

U = [ρ, ρu, ρv, ρE]ᵀ
F = [ρu, ρu² + p, ρuv, u(ρE + p)]ᵀ
G = [ρv, ρuv, ρv² + p, v(ρE + p)]ᵀ

p = (γ - 1)·(ρE - 0.5·ρ·(u² + v²))
```

Expected result: Cp distribution with a shock on the upper surface
somewhere around x/c ≈ 0.55–0.60. For a pure inviscid case the shock
will be stronger and slightly further aft than the viscous experiment
(AGARD AR-138), which is fine — do not tune the solver to match the
viscous data.

---

## 2. Project structure

Create this layout. Keep it flat and readable, no frameworks.

```
euler2d/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── Airfoil.hpp      // RAE 2822 coordinate table + interpolation
│   ├── Mesh.hpp         // C-grid struct + generator + I/O
│   ├── Solver.hpp       // Main solver class
│   ├── Flux.hpp         // Roe flux + MUSCL reconstruction
│   ├── BC.hpp           // Boundary conditions (ghost cells)
│   └── IO.hpp           // Tecplot/VTK writers, Cp extraction
├── src/
│   ├── main.cpp         // Parses config, runs the pipeline
│   ├── Mesh.cpp
│   ├── Solver.cpp
│   ├── Flux.cpp
│   ├── BC.cpp
│   └── IO.cpp
├── config/
│   └── rae2822.ini      // Runtime parameters
└── output/              // Generated at runtime (gitignored)
    ├── mesh.dat         // Tecplot ASCII — mesh only (for debug)
    ├── field.vts        // VTK StructuredGrid — full flow field for ParaView
    ├── surface.vtp      // VTK PolyData — airfoil surface with Cp
    ├── cp.dat           // ASCII (x/c, Cp) — for plotting
    └── convergence.dat  // ASCII (iter, residual)
```

Config file format (simple key=value):

```ini
# rae2822.ini
mach        = 0.729
aoa_deg     = 2.31
gamma       = 1.4

ni_wrap     = 257      # points around airfoil + wake
nj_normal   = 65       # points in normal direction
n_airfoil   = 129      # points on airfoil (half+half-1)
n_wake      = 64       # points in each wake branch
r_far       = 25.0     # farfield radius in chords

cfl         = 1.5
max_iter    = 20000
residual_drop = 1.0e-5   # stop when ||ρ residual|| drops this much
scheme_order = 2         # 1 = first order, 2 = MUSCL
output_interval = 500    # write field.vts / surface.vtp every N iterations
```

---

## 3. Numerical methods — exact specification

### 3.1 Grid: C-topology

- Cell-centered finite volume. Store cell centers, cell areas, and
  face normals + lengths.
- Index convention: `i ∈ [0, ni)` wraps from lower wake outflow →
  lower airfoil surface → LE → upper airfoil surface → upper wake
  outflow. `j ∈ [0, nj)` goes from wall/wake (j=0) to farfield (j=nj-1).
- Generate the airfoil curve via cosine clustering:
  `x_k = 0.5·(1 - cos(π·k/(N-1)))` for k = 0..N-1.
- Use linear interpolation into the RAE 2822 coordinate table (hardcode
  it in `Airfoil.hpp` — do NOT download anything).
- Outer boundary: C-shape. Two straight segments at `x = r_far` extending
  from `y = ±r_far` to `y = 0`, joined by a half-circle of radius `r_far`
  centered at the origin on the `x < 0` side.
- Fill the domain interior by transfinite interpolation (TFI), linear
  blending between inner (airfoil+wake) and outer (C-shape) curves,
  with tanh clustering in j-direction for near-wall resolution.

### 3.2 Finite-volume discretization

Cell-centered. For each cell `(i,j)` with area `A_ij` and four faces:

```
dU_ij/dt = -(1/A_ij) · Σ_faces (F_n · Δs)
```

Store conserved variables `U[i][j][4]` at cell centers (ghost cells
included).

### 3.3 Roe flux with Harten entropy fix

For each face, given left/right states `U_L, U_R`:

```
F_Roe = 0.5·(F(U_L) + F(U_R)) - 0.5·|A_roe|·(U_R - U_L)
```

where `A_roe` is evaluated at Roe-averaged state:

```
ρ̂  = sqrt(ρ_L · ρ_R)
û   = (sqrt(ρ_L)·u_L + sqrt(ρ_R)·u_R) / (sqrt(ρ_L) + sqrt(ρ_R))
v̂   = (sqrt(ρ_L)·v_L + sqrt(ρ_R)·v_R) / (sqrt(ρ_L) + sqrt(ρ_R))
Ĥ   = (sqrt(ρ_L)·H_L + sqrt(ρ_R)·H_R) / (sqrt(ρ_L) + sqrt(ρ_R))
â   = sqrt((γ-1)·(Ĥ - 0.5·(û² + v̂²)))
```

Eigenvalues: `λ1 = q̂ - â, λ2 = λ3 = q̂, λ4 = q̂ + â`
where `q̂ = û·nx + v̂·ny` is the face-normal velocity.

**Harten entropy fix** (critical — otherwise expansion shocks will
appear and the solver can diverge or produce garbage at transonic
conditions):

```
if |λ_k| < δ:   |λ_k| ← (λ_k² + δ²) / (2·δ)
```

with `δ = 0.1·â` typically. Apply only to acoustic eigenvalues
(λ1, λ4). Leaving the contact eigenvalues alone is fine.

Rotate to face-normal frame, compute the flux, rotate back. Write this
in `Flux.cpp` — this is the single most error-prone piece. Test it
alone on a supersonic contact discontinuity before plugging into the
solver.

### 3.4 MUSCL reconstruction (2nd order)

For face between cells `i` and `i+1`:

```
U_L = U_i   + 0.5·φ(r_i)  ·(U_i   - U_{i-1})
U_R = U_{i+1} - 0.5·φ(r_{i+1})·(U_{i+2} - U_{i+1})
```

Use **van Albada limiter**:
```
φ(r) = (r² + r) / (r² + 1),   r = (U_{i+1} - U_i) / (U_i - U_{i-1})
```

Reconstruct on **primitive variables** (ρ, u, v, p), not conservative
— this is more robust near shocks.

### 3.5 Boundary conditions (ghost-cell approach)

Add one layer of ghost cells around the domain (j = -1 and j = nj for
far/wall, plus wrap handling for i).

**Solid wall (j = 0, airfoil surface segment):** Reflected-velocity
ghost cell.
```
ρ_ghost = ρ_interior
u_ghost = u_interior - 2·(u_interior · n̂)·n̂_x
v_ghost = v_interior - 2·(u_interior · n̂)·n̂_y
p_ghost = p_interior
```
This makes the normal velocity at the face zero and leaves tangential
untouched.

**Farfield (j = nj-1):** Characteristic BC using Riemann invariants.
Compute the normal velocity `q_n`, local sound speed `a`, and Mach
number. Decide based on sign of `q_n` and whether `|M_n| < 1`:
- Supersonic inflow: copy all 4 freestream values
- Supersonic outflow: extrapolate all 4 from interior
- Subsonic inflow: 3 from freestream, 1 from interior (outgoing
  characteristic `R⁺ = q_n + 2a/(γ-1)`)
- Subsonic outflow: 1 from freestream (incoming `R⁻ = q_n - 2a/(γ-1)`),
  3 from interior (entropy, tangential velocity, `R⁺`)

Most cells on the farfield of this case will be **subsonic** (freestream
M = 0.729 < 1). Handle it correctly — this is the second most
error-prone thing after the Roe flux.

**Wake cut (j = 0, i ∈ [0, i_TE_lower) ∪ (i_TE_upper, ni)):**
The two wake branches of the C-grid are geometrically coincident.
Treat them as an **internal interface**: the ghost cell for the lower
wake point at index `i` is the interior cell across the cut from the
upper wake. Concretely, for a lower-wake cell at `(i, 0)` its ghost
neighbor across the j=0 face is the upper-wake interior cell at
`(i_mirror, 0)` where `i_mirror = ni - 1 - i` (or whatever matches
your indexing — work out and document the exact formula).

### 3.6 Time integration

Explicit 3-stage SSP Runge-Kutta:
```
U(1) = U^n + Δt·R(U^n)
U(2) = 0.75·U^n + 0.25·(U(1) + Δt·R(U(1)))
U^{n+1} = (1/3)·U^n + (2/3)·(U(2) + Δt·R(U(2)))
```

Use **local time stepping**: for each cell,
`Δt_ij = CFL · A_ij / (Σ |λ_max,face| · Δs_face)`.

### 3.7 Convergence monitoring

L2 norm of density residual:
```
res = sqrt( Σ (R_ρ)² / N_cells )
```
Normalize by the residual at iteration 1. Stop when it drops by
`residual_drop` (5 orders typically) or `max_iter` reached.

### 3.8 Output — field export to ParaView

Export must be openable directly in ParaView without any conversion.
Use **VTK XML formats** (ASCII is fine, no need for binary/base64 —
this is a small 2D grid, human-readable is more useful for debugging).

**`field.vts` — StructuredGrid** of the full flow field. This is the
primary output the user will inspect after each run.

Node layout: write the grid nodes (cell *corners*, not centers) of the
C-grid. Cell-centered data is averaged to nodes for export, OR write
the cell-centered values as CellData (preferred — avoids averaging
artifacts near the shock). Use `<CellData>` for:

| Field name    | Unit / definition                                    |
|---------------|------------------------------------------------------|
| `density`     | ρ (non-dimensional, ρ∞ = 1)                          |
| `pressure`    | p                                                    |
| `velocity`    | vector (u, v, 0)  — write as 3-component for VTK     |
| `mach`        | \|u\| / a, where a = sqrt(γ·p/ρ)                     |
| `cp`          | (p - p∞) / (0.5·ρ∞·\|u∞\|²)                          |
| `temperature` | T = p / (ρ·R), with R = 1 in non-dim (T = γ·p/ρ)     |
| `entropy`     | s = ln(p/ρ^γ) — useful for checking entropy fix      |

VTK structured grid format (`.vts`) — minimal skeleton:

```xml
<?xml version="1.0"?>
<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">
  <StructuredGrid WholeExtent="0 NI 0 NJ 0 0">
    <Piece Extent="0 NI 0 NJ 0 0">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          x0 y0 0  x1 y1 0  ...
        </DataArray>
      </Points>
      <CellData Scalars="density" Vectors="velocity">
        <DataArray type="Float64" Name="density" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="pressure" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="mach" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="cp" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="temperature" format="ascii"> ... </DataArray>
        <DataArray type="Float64" Name="entropy" format="ascii"> ... </DataArray>
      </CellData>
    </Piece>
  </StructuredGrid>
</VTKFile>
```

Notes on VTK conventions:
- `Extent` uses **node counts minus 1** (i.e. number of cells in each
  direction). For a grid with NI nodes in i and NJ nodes in j, extent
  is `"0 NI-1 0 NJ-1 0 0"`.
- Point ordering: fastest index is i, then j, then k. So loop
  `for j: for i: write(x[i,j], y[i,j], 0)`.
- CellData ordering matches point ordering but over cells — same loop
  pattern, one less entry per direction.
- Do **not** include ghost cells in the export. Export only real cells
  `i ∈ [0, ni-1)`, `j ∈ [0, nj-1)` (where ni, nj are node counts).

**`surface.vtp` — PolyData** of the airfoil surface (a 1D polyline
in 2D space), with per-point `cp`, `pressure`, `mach_surface`.
This lets you color the airfoil surface by Cp directly in ParaView
alongside the field. Extract along `j = 0` between `i_TEl` and
`i_TEu`.

**`cp.dat`** — plain two-column ASCII (x/c, Cp) along the airfoil
surface, sorted by i-index. Easy to plot in gnuplot/Python.

**`convergence.dat`** — three columns: iteration, absolute L2 residual
of density, residual normalized by iter-1 value.

**`mesh.dat`** — Tecplot ASCII of the mesh only (nodes, no field).
This is written in Stage 1 before any solver exists, purely for visual
verification of the grid generator.

**Output frequency:** write `convergence.dat` every iteration (just
append). Write `field.vts`, `surface.vtp`, `cp.dat` every
`output_interval` iterations (add this to the INI, default 500) and
once at the end.

---

## 4. Acceptance criteria

The build succeeds only when **all** of these hold:

1. `cmake --build build` compiles without warnings using `-Wall -Wextra`.
2. Running `./euler2d config/rae2822.ini` produces, in `output/`:
   - `mesh.dat` — Tecplot ASCII of the mesh, opens cleanly in ParaView
     (no negative cell areas, no crossed cells; verify by eye).
   - `field.vts` — VTK StructuredGrid with CellData: density, pressure,
     velocity, mach, cp, temperature, entropy. Must open in ParaView
     in one click and render without warnings.
   - `surface.vtp` — VTK PolyData of the airfoil surface with Cp.
   - `cp.dat` — two columns (x/c, Cp) along the airfoil surface only
     (j=0 cells between iTEl and iTEu).
   - `convergence.dat` — (iteration, residual, normalized residual).
3. Density residual drops at least 4 orders of magnitude.
4. Cp plot shows a recognizable transonic shock on the upper surface
   between x/c = 0.5 and x/c = 0.7.
5. Lift coefficient CL is in the rough range [0.7, 0.95] (the viscous
   experiment gives ~0.743; inviscid will be somewhat higher).
6. No NaNs anywhere in the field.
7. In ParaView, applying a "Contour" filter on `mach` with a threshold
   of 1.0 must produce a visible sonic line wrapping the upper surface
   and terminating at the shock.

---

## 5. Execution plan — do this strictly in order

Do **not** write everything up front. Implement one stage, build it,
verify the stated checkpoint, commit, then move on. If a checkpoint
fails, fix before proceeding. Ask me before making any substantive
architectural deviation from this plan.

### Stage 1 — Mesh only
- Implement `Airfoil.hpp` and `Mesh.hpp/cpp`.
- Write a `main.cpp` that only generates the mesh and dumps Tecplot.
- **Checkpoint:** open `mesh.dat` in ParaView. Airfoil shape is correct,
  C-grid has no self-intersections, cells near the wall are thin and
  smooth, wake cut is properly collapsed on y=0 behind the TE.
  STOP and show me the visualization before continuing.

### Stage 2 — Data structures and BCs (no flux yet)
- Cell-centered storage, ghost cells, metrics (areas, face normals).
- Initialize the whole field to freestream.
- Apply BCs once and dump the field.
- **Checkpoint:** ghost cells have sensible values; wall ghost cells
  show reflected normal velocity; wake-cut ghost cells correctly
  mirror their partner. Print a few sample ghost values and let me
  verify.

### Stage 3 — First-order Roe solver
- Implement Roe flux with entropy fix in `Flux.cpp`.
- Unit-test it: on a uniform freestream state, the flux should reduce
  to the exact freestream flux (zero dissipation contribution).
  On a stationary contact, mass/momentum/energy should be consistent.
- Wire into the solver with piecewise-constant reconstruction
  (no MUSCL yet).
- Implement the VTK writer (`field.vts`) now — we want to visualize
  from the first solver run, not at the end.
- Run to steady state.
- **Checkpoint:** residual drops ≥ 3 orders, shock appears but is
  smeared over ~10 cells. `field.vts` opens in ParaView and shows
  a recognizable transonic flow pattern: acceleration over the upper
  surface, supersonic pocket with M > 1, shock. Show me a screenshot
  of Mach contours from ParaView.

### Stage 4 — Second-order MUSCL
- Add primitive-variable MUSCL + van Albada.
- **Checkpoint:** residual drops ≥ 4–5 orders, shock sharpens to
  ~2–3 cells, no oscillations upstream of shock.

### Stage 5 — Output and comparison
- Add `surface.vtp` writer for the airfoil surface.
- Produce final Cp plot (gnuplot or similar from `cp.dat`).
- Show me: (a) Mach contour screenshot from ParaView, (b) Cp plot,
  (c) convergence history.

---

## 6. Rules of engagement

- **Language:** C++17. No external libraries except the standard
  library. No Eigen, no Boost, no nothing. This stays self-contained.
- **Build:** single CMakeLists.txt, target `euler2d`. Default to
  Release with `-O2 -march=native`.
- **Style:** prefer clarity over cleverness. Use `double` everywhere.
  Put the physics-heavy functions (Roe flux, entropy fix, char BC) in
  clearly-commented blocks with references to the equations above.
- **No silent decisions.** If anything in this spec is ambiguous or if
  you hit a choice point not covered (e.g., exact form of entropy fix
  constant, CFL for Forward Euler vs RK3, how to handle sharp-TE
  singularity), ask before committing to a choice.
- **Commit after each stage** with a message like
  `stage 1: C-grid generator + Tecplot export`.

Begin with Stage 1. Ask me anything unclear before writing code.
