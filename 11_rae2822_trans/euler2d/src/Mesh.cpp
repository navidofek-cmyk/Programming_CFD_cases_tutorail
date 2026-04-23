#include "Mesh.hpp"
#include "Airfoil.hpp"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <cassert>

static constexpr double PI = 3.14159265358979323846;

// ---- helpers ---------------------------------------------------------------

static double tanh_stretch(double t, double beta) {
    // Cluster near t=0 (wall): s=0 at t=0, s=1 at t=1, dense near wall.
    return 1.0 - std::tanh(beta * (1.0 - t)) / std::tanh(beta);
}

// ---- inner boundary (j=0): airfoil + wake lines ----------------------------
//
// C-grid i-index layout (ni_wrap nodes total):
//
//   [0 .. n_wake-1]                    : lower wake (x goes from r_far → TE)
//   [n_wake .. n_wake + n_airfoil - 1] : airfoil (lower TE → LE → upper TE)
//                                         sub-index 0 = lower TE, mid = LE, end = upper TE
//   [n_wake + n_airfoil - 1 .. ni-1]   : upper wake (TE → r_far)
//
// So:
//   i_TEl = n_wake - 1      (lower TE, also j=0 end of lower-wake branch)
//           actually node n_wake is the first airfoil node = lower TE
//   Let's define:
//   i_TEl = n_wake           (lower trailing edge)
//   i_LE  = n_wake + (n_airfoil-1)/2
//   i_TEu = n_wake + n_airfoil - 1  (upper trailing edge)
//
// Wake nodes i in [0, n_wake-1] and [i_TEu+1, ni-1] lie on y=0, x > 1.

static void build_inner_boundary(Mesh& m, const MeshConfig& cfg) {
    int nw = cfg.n_wake;
    int na = cfg.n_airfoil;  // must be odd so LE lands on a node
    assert((na % 2) == 1);  // n_airfoil must be odd

    m.i_TEl = nw;
    m.i_LE  = nw + (na - 1) / 2;
    m.i_TEu = nw + na - 1;

    double r_far = cfg.r_far;

    // --- lower wake: i in [0, nw-1], y=0, x from r_far down to 1.0 ---
    // cosine clustering so that spacing near TE matches airfoil TE spacing
    for (int i = 0; i < nw; ++i) {
        double t = (double)i / (nw - 1);              // 0 at i=0 (far), 1 at i=nw-1 (TE)
        double s = 0.5 * (1.0 - std::cos(PI * t));    // cosine: fine near both ends
        double xw = r_far + s * (1.0 - r_far);        // r_far → 1
        m.node_x(i, 0) = xw;
        m.node_y(i, 0) = 0.0;
    }

    // --- airfoil: i in [nw, nw+na-1] ---
    // Half-half split: [nw .. i_LE] = lower surface (TE→LE),
    //                  [i_LE .. i_TEu] = upper surface (LE→TE)
    int n_half = (na - 1) / 2;  // nodes on each half, not counting LE

    // Lower surface: k=0 → TE (x=1), k=n_half → LE (x=0)
    for (int k = 0; k <= n_half; ++k) {
        double t  = (double)k / n_half;
        double xc = 0.5 * (1.0 - std::cos(PI * t));  // 0 at LE, 1 at TE
        double xl = 1.0 - xc;                          // TE=1→LE=0 (going LE side)
        // k=0 is TE (xl=1), k=n_half is LE (xl=0)
        int gi = m.i_TEl + k;  // global i index
        m.node_x(gi, 0) = xl;
        m.node_y(gi, 0) = airfoil_lower(xl);
    }

    // Upper surface: k=0 → LE (x=0), k=n_half → TE (x=1)
    for (int k = 0; k <= n_half; ++k) {
        double t  = (double)k / n_half;
        double xc = 0.5 * (1.0 - std::cos(PI * t));  // 0→1
        int gi = m.i_LE + k;
        m.node_x(gi, 0) = xc;
        m.node_y(gi, 0) = airfoil_upper(xc);
    }
    // LE node set twice (once from lower, once from upper) — should match exactly

    // --- upper wake: i in [i_TEu, ni-1], y=0, x from 1 to r_far ---
    for (int k = 0; k <= nw - 1; ++k) {
        double t = (double)k / (nw - 1);
        double s = 0.5 * (1.0 - std::cos(PI * t));
        double xw = 1.0 + s * (r_far - 1.0);
        int gi = m.i_TEu + k;
        m.node_x(gi, 0) = xw;
        m.node_y(gi, 0) = 0.0;
    }
}

// ---- outer boundary (j=nj-1): C-shape at radius r_far ---------------------
//
// Parametrize with the same i-index as inner boundary.
// The outer boundary is:
//   lower straight segment: x = r_far, y from 0 to -r_far  (lower wake side)
//   semicircle: angle from -pi/2 to +pi/2 (going CCW) around origin
//   upper straight segment: x = r_far, y from r_far to 0
//
// We match the i-distribution of the inner boundary in the same proportions.

static void build_outer_boundary(Mesh& m, const MeshConfig& cfg) {
    int nw = cfg.n_wake;
    int na = cfg.n_airfoil;
    double r = cfg.r_far;
    int nj = m.nj;

    // Outer C-boundary — three pieces, junction nodes shared:
    //
    //  Lower vertical:  i=[0, nw-1],       x=r, y: 0 → -r
    //  Arc:             i=[nw-1, nw+na-1],  (r,-r) → (-r,0) → (r,r)
    //  Upper vertical:  i=[nw+na-1, ni-1],  x=r, y: r → 0
    //
    // The arc lies on a circle centered at (r/4, 0), radius 5r/4.
    // That circle passes through (r,-r), (-r,0), and (r,r), giving a
    // smooth D-shaped farfield with no diagonal "triangle" in the wake.

    // Lower vertical: x=r, y from 0 to -r
    for (int i = 0; i < nw; ++i) {
        double t = (double)i / (nw - 1);
        m.node_x(i, nj - 1) = r;
        m.node_y(i, nj - 1) = -r * t;
    }
    // i=nw-1: (r, -r)

    // Arc on circle (cx=r/4, cy=0, R=5r/4):
    // theta_start = atan2(-r, r-r/4) = atan2(-4,3), going CW through -pi to theta_end
    double cx = r / 4.0;
    double R  = 5.0 * r / 4.0;
    double theta_start = std::atan2(-r, r - cx);           // ≈ -0.9273 rad
    double theta_end   = std::atan2( r, r - cx);           // ≈ +0.9273 rad
    double arc_span    = 2.0 * PI - (theta_end - theta_start);  // long way around

    for (int k = 0; k <= na; ++k) {
        double t = (double)k / na;
        double theta = theta_start - t * arc_span;          // CW → goes via x<0
        int gi = (nw - 1) + k;
        m.node_x(gi, nj - 1) = cx + R * std::cos(theta);
        m.node_y(gi, nj - 1) =      R * std::sin(theta);
    }
    // gi=nw-1 (k=0): (r,-r)  — same as lower vertical end ✓
    // gi=nw+na-1 (k=na): (r,r)

    // Upper vertical: x=r, y from r to 0
    for (int k = 0; k < nw; ++k) {
        double t = (double)k / (nw - 1);
        int gi = (nw + na - 1) + k;
        m.node_x(gi, nj - 1) = r;
        m.node_y(gi, nj - 1) = r * (1.0 - t);
    }
    // gi=nw+na-1 (k=0): (r,r)  — same as arc end ✓
    // gi=ni-1 (k=nw-1): (r,0)
}

// ---- TFI fill interior ------------------------------------------------------
//
// Linear TFI between inner (j=0) and outer (j=nj-1) curves,
// with tanh clustering in j to pack cells near the wall.

static void tfi_fill(Mesh& m, const MeshConfig& /*cfg*/) {
    int ni = m.ni;
    int nj = m.nj;
    double beta = 3.5;  // tanh stretching, clusters near wall (j=0)

    for (int i = 0; i < ni; ++i) {
        double x0 = m.node_x(i, 0);
        double y0 = m.node_y(i, 0);
        double x1 = m.node_x(i, nj - 1);
        double y1 = m.node_y(i, nj - 1);

        for (int j = 0; j < nj; ++j) {
            double t = (double)j / (nj - 1);
            double s = tanh_stretch(t, beta);
            m.node_x(i, j) = x0 + s * (x1 - x0);
            m.node_y(i, j) = y0 + s * (y1 - y0);
        }
    }
}

// ---- compute_metrics --------------------------------------------------------

void Mesh::compute_metrics() {
    int nci = nc_i();
    int ncj = nc_j();

    xc.resize(nci * ncj);
    yc.resize(nci * ncj);
    area.resize(nci * ncj);

    // Cell centers (average of 4 nodes)
    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
            double cx = 0.25 * (node_x(i,j) + node_x(i+1,j) + node_x(i+1,j+1) + node_x(i,j+1));
            double cy = 0.25 * (node_y(i,j) + node_y(i+1,j) + node_y(i+1,j+1) + node_y(i,j+1));
            xc[j * nci + i] = cx;
            yc[j * nci + i] = cy;

            // Cell area via cross product of diagonals (shoelace for quad)
            // Nodes: (i,j)=A, (i+1,j)=B, (i+1,j+1)=C, (i,j+1)=D
            double ax = node_x(i,j),   ay = node_y(i,j);
            double bx = node_x(i+1,j), by = node_y(i+1,j);
            double cx2= node_x(i+1,j+1), cy2= node_y(i+1,j+1);
            double dx = node_x(i,j+1), dy = node_y(i,j+1);
            // Shoelace: A=0.5*|AC×BD| for quad ABCD with diagonals AC and BD
            double ac_x = cx2 - ax, ac_y = cy2 - ay;
            double bd_x = dx  - bx, bd_y = dy  - by;
            double A = 0.5 * std::abs(ac_x * bd_y - ac_y * bd_x);
            area[j * nci + i] = A;
        }
    }

    // i-faces (between column i and i+1, for j in [0,ncj-1])
    // face runs from node(i+1,j) to node(i+1,j+1) — wait, actually:
    // i-face at (i+1/2, j) separates cell (i,j) and (i+1,j).
    // Face nodes: (i+1, j) and (i+1, j+1).
    // Normal points in +i direction (from cell i to i+1).
    // n = (dy, -dx)/len where (dx,dy) = face vector along j.
    ifnx.resize(ni * ncj);
    ifny.resize(ni * ncj);
    iflen.resize(ni * ncj);

    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i < ni; ++i) {
            // Face between node (i,j) and (i,j+1)
            double dx = node_x(i, j+1) - node_x(i, j);
            double dy = node_y(i, j+1) - node_y(i, j);
            double len = std::sqrt(dx*dx + dy*dy);
            // Outward normal from cell (i-1,j) to cell (i,j): (+dy, -dx)/len
            ifnx[j * ni + i] =  dy / len;
            ifny[j * ni + i] = -dx / len;
            iflen[j * ni + i] = len;
        }
    }

    // j-faces (between row j and j+1)
    // Face at (i, j+1/2) separates cell (i,j) and (i,j+1).
    // Face nodes: (i,j+1) and (i+1,j+1).
    // Normal in +j direction.
    jfnx.resize(nci * nj);
    jfny.resize(nci * nj);
    jflen.resize(nci * nj);

    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < nci; ++i) {
            double dx = node_x(i+1, j) - node_x(i, j);
            double dy = node_y(i+1, j) - node_y(i, j);
            double len = std::sqrt(dx*dx + dy*dy);
            // Outward normal from cell (i,j-1) to cell (i,j): (-dy, dx)/len
            jfnx[j * nci + i] = -dy / len;
            jfny[j * nci + i] =  dx / len;
            jflen[j * nci + i] = len;
        }
    }
}

// ---- generate ---------------------------------------------------------------

void Mesh::generate(const MeshConfig& cfg) {
    if ((cfg.n_airfoil % 2) == 0)
        throw std::invalid_argument("n_airfoil must be odd");

    ni = cfg.ni_wrap;
    nj = cfg.nj_normal;

    // ni_wrap = 2*n_wake + n_airfoil - 1
    // (the TE nodes at both ends of the airfoil are shared with the wake)
    int expected_ni = 2 * cfg.n_wake + cfg.n_airfoil - 1;
    if (ni != expected_ni)
        throw std::invalid_argument(
            "ni_wrap must equal 2*n_wake + n_airfoil - 1");

    x.assign(ni * nj, 0.0);
    y.assign(ni * nj, 0.0);

    build_inner_boundary(*this, cfg);
    build_outer_boundary(*this, cfg);
    tfi_fill(*this, cfg);
    compute_metrics();
}

// ---- write_tecplot ----------------------------------------------------------

void Mesh::write_tecplot(const std::string& filename) const {
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot open " + filename);

    f << "TITLE = \"RAE2822 C-grid\"\n";
    f << "VARIABLES = \"X\" \"Y\"\n";
    f << "ZONE T=\"Mesh\" I=" << ni << " J=" << nj
      << " DATAPACKING=POINT\n";

    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            f << node_x(i, j) << " " << node_y(i, j) << "\n";
        }
    }
}
