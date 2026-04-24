#include "BC.hpp"
#include "Mesh.hpp"
#include <cmath>

// ---- i-exit ghosts (i=-1 and i=nci): zeroth-order extrapolation -----------
// These are the far-wake outflow boundaries at x=r_far, y≈0.
void bc_i_exit(Solver& s) {
    int nci = s.nci;
    int ncj = s.ncj;
    for (int j = 0; j < ncj; ++j) {
        for (int var = 0; var < 4; ++var) {
            s.U(var, -1,  j) = s.U(var, 0,      j);
            s.U(var, nci, j) = s.U(var, nci - 1, j);
        }
    }
}

// ---- wall BC (j=-1 ghost) ---------------------------------------------------
//
// Inviscid wall: no normal velocity through wall.
// Reflect the interior cell (j=0) across the wall face:
//   u_ghost = u_real - 2*(u·n)*n
// where n is the inward wall normal (jfn at j=0 row points outward in +j,
// so inward wall normal = -jfn).
// Density and pressure are mirrored unchanged.

void bc_wall(Solver& s) {
    const Mesh& m = *s.pmesh;
    int nci = s.nci;

    // j-face normals at j=0 (separates wall ghost j=-1 from real cell j=0)
    // jfnx[0*nci + i], jfny[0*nci + i] point in +j direction (away from wall)
    for (int i = 0; i < nci; ++i) {
        double nx =  m.jfnx[i];   // outward from ghost (= inward wall normal in +j)
        double ny =  m.jfny[i];

        double r  = s.rho  (i, 0);
        double u  = s.vel_u(i, 0);
        double v  = s.vel_v(i, 0);
        double p  = s.press(i, 0);

        // normal component of velocity
        double un = u * nx + v * ny;

        // reflected velocity: subtract twice normal component
        double u_g = u - 2.0 * un * nx;
        double v_g = v - 2.0 * un * ny;

        s.U(0, i, -1) = r;
        s.U(1, i, -1) = r * u_g;
        s.U(2, i, -1) = r * v_g;
        // Total energy from mirrored velocity + same pressure
        s.U(3, i, -1) = p / (s.gamma - 1.0) + 0.5 * r * (u_g*u_g + v_g*v_g);
    }
}

// ---- farfield BC (j=ncj ghost) -----------------------------------------------
//
// Riemann characteristic BC for subsonic farfield.
// One Riemann invariant comes from inside, one from freestream.
// R+ = u_n + 2c/(γ-1)  propagates outward (from domain to farfield)
// R- = u_n - 2c/(γ-1)  propagates inward  (from farfield to domain)
//
// For outflow: u_n > 0 (normal pointing outward = +j)
//   R+ from interior, R- from freestream → u_n_g, c_g
// For inflow: u_n < 0
//   R+ from freestream, R- from interior → same formula with swapped source

void bc_farfield(Solver& s) {
    const Mesh& m = *s.pmesh;
    int nci = s.nci;
    int ncj = s.ncj;

    double g   = s.gamma;
    double gm1 = g - 1.0;

    // Freestream
    double r_inf = s.rho_inf;
    double u_inf = s.u_inf;
    double v_inf = s.v_inf;
    double p_inf = s.p_inf;
    double c_inf = std::sqrt(g * p_inf / r_inf);

    // j-face at j=ncj: jfnx[ncj*nci + i], jfny[ncj*nci + i]
    // These normals point in +j direction (outward from domain)
    for (int i = 0; i < nci; ++i) {
        double nx = m.jfnx[ncj * nci + i];
        double ny = m.jfny[ncj * nci + i];

        // Interior cell is j=ncj-1
        double r_i = s.rho  (i, ncj - 1);
        double u_i = s.vel_u(i, ncj - 1);
        double v_i = s.vel_v(i, ncj - 1);
        double p_i = s.press(i, ncj - 1);
        double c_i = std::sqrt(g * p_i / r_i);

        double un_i   = u_i * nx + v_i * ny;
        double un_inf = u_inf * nx + v_inf * ny;

        double Rp, Rm;
        if (un_i >= 0.0) {
            // outflow: R+ from inside, R- from freestream
            Rp = un_i   + 2.0 * c_i   / gm1;
            Rm = un_inf - 2.0 * c_inf / gm1;
        } else {
            // inflow: both invariants from freestream
            Rp = un_inf + 2.0 * c_inf / gm1;
            Rm = un_inf - 2.0 * c_inf / gm1;
        }

        double un_g = 0.5 * (Rp + Rm);
        double c_g  = 0.25 * gm1 * (Rp - Rm);
        if (c_g < 1e-10) c_g = 1e-10;

        double s_i = p_i / std::pow(r_i, g);  // isentropic entropy
        double r_g = std::pow(c_g * c_g / (g * s_i), 1.0 / gm1);
        double p_g = r_g * c_g * c_g / g;

        // Tangential velocity from upstream (freestream for inflow, interior for outflow)
        double ut_x, ut_y;
        if (un_i >= 0.0) {
            // outflow: tangential from interior
            ut_x = u_i - un_i * nx;
            ut_y = v_i - un_i * ny;
        } else {
            // inflow: tangential from freestream
            ut_x = u_inf - un_inf * nx;
            ut_y = v_inf - un_inf * ny;
        }

        double u_g = ut_x + un_g * nx;
        double v_g = ut_y + un_g * ny;

        s.U(0, i, ncj) = r_g;
        s.U(1, i, ncj) = r_g * u_g;
        s.U(2, i, ncj) = r_g * v_g;
        s.U(3, i, ncj) = p_g / gm1 + 0.5 * r_g * (u_g*u_g + v_g*v_g);
    }
}

// ---- wake-cut BC (j=-1 ghost for wake cells) ---------------------------------
//
// The C-grid has a wake cut along y=0 for x > 1 (trailing edge wake).
// The lower branch (i in [0, i_TEl-1]) and upper branch (i in [i_TEu, nci-1])
// are geometrically coincident but topologically distinct.
// Ghost cell j=-1 of cell i in the lower wake is the mirror of real cell
// j=0 of cell (nci-1-i) in the upper wake, and vice-versa.
//
// C-grid cell indexing (nci=255):
//   i in [0  .. i_TEl-1] = [0..63]   : lower-wake cells
//   i in [i_TEl..i_TEu-1]= [64..191] : airfoil cells (wall BC applied above)
//   i in [i_TEu..nci-1]  = [192..254]: upper-wake cells
//
// Mirror map: i_mirror(i) = nci - 1 - i
//   i=0   ↔ i=254  (lower-wake far ↔ upper-wake far)
//   i=63  ↔ i=191  (lower-wake TE  ↔ upper-wake TE start)
//   i=192 ↔ i=63   (upper-wake TE  ↔ lower-wake TE end)

void bc_wake_cut(Solver& s) {
    const Mesh& m = *s.pmesh;
    int nci = s.nci;
    int i_TEl = m.i_TEl;    // = n_wake (first airfoil node index = 64)
    int i_TEu = m.i_TEu;    // = n_wake + n_airfoil - 1 (= 192)

    // lower wake: i in [0, i_TEl-1]
    for (int i = 0; i < i_TEl; ++i) {
        int im = nci - 1 - i;  // mirror index in upper wake
        for (int var = 0; var < 4; ++var)
            s.U(var, i, -1) = s.U(var, im, 0);
    }

    // upper wake: i in [i_TEu, nci-1]
    for (int i = i_TEu; i < nci; ++i) {
        int im = nci - 1 - i;  // mirror index in lower wake
        for (int var = 0; var < 4; ++var)
            s.U(var, i, -1) = s.U(var, im, 0);
    }
}
