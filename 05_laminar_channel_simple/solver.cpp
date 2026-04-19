#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary.hpp"

// ---------------------------------------------------------------------------
// Momentum predictor — u component
// Interior u faces: i in [1, Nx-1], j in [0, Ny-1]
//
// No-slip at top/bottom walls is handled via ghost cell:
//   u_ghost_south = -u(i, 0)   (enforces u=0 at y=0)
//   u_ghost_north = -u(i,Ny-1) (enforces u=0 at y=Ly)
// This adds the wall contribution to aP (see comment in loop).
// ---------------------------------------------------------------------------
static void u_momentum(const Mesh& m, const SolverSettings& s,
                        const Field& u_old, const Field& v_old,
                        const Field& p,
                        Field& u_star, Field& aP_u) {
    const double dx = m.dx, dy = m.dy;
    const double De = s.nu * dy / dx;
    const double Dn = s.nu * dx / dy;

    for (int j = 0; j < m.Ny; ++j) {
        for (int i = 1; i < m.Nx; ++i) {

            // Face-centred advecting velocities
            const double u_e = 0.5 * (u_old(i, j) + u_old(i + 1, j));
            const double u_w = 0.5 * (u_old(i - 1, j) + u_old(i, j));
            // v at north/south faces of the u-CV; wall faces have v=0
            const double v_n = (j < m.Ny - 1) ? 0.5 * (v_old(i - 1, j + 1) + v_old(i, j + 1)) : 0.0;
            const double v_s = (j > 0)         ? 0.5 * (v_old(i - 1, j)     + v_old(i, j))     : 0.0;

            const double Fe = s.rho * u_e * dy;
            const double Fw = s.rho * u_w * dy;
            const double Fn = s.rho * v_n * dx;
            const double Fs = s.rho * v_s * dx;

            const double aE  = De + std::max(-Fe, 0.0);
            const double aW  = De + std::max( Fw, 0.0);
            const double aN  = Dn + std::max(-Fn, 0.0);
            const double aS  = Dn + std::max( Fs, 0.0);
            const double aP0 = s.rho * dx * dy / s.pseudo_dt;
            const double aP  = aE + aW + aN + aS + aP0;

            // Pressure source (exact on staggered grid — no averaging!)
            const double b = -(p(i, j) - p(i - 1, j)) * dy;

            double rhs    = aP0 * u_old(i, j) + b;
            rhs += aE * u_old(i + 1, j);  // outlet BC already set in u_old(Nx,j)
            rhs += aW * u_old(i - 1, j);  // inlet BC: u_old(0,j) = u_in

            double aP_eff = aP;
            // North neighbour (or wall ghost = -u_P)
            if (j < m.Ny - 1) {
                rhs += aN * u_old(i, j + 1);
            } else {
                aP_eff += aN;  // ghost: -u_P → contribution moves to diagonal
            }
            // South neighbour (or wall ghost = -u_P)
            if (j > 0) {
                rhs += aS * u_old(i, j - 1);
            } else {
                aP_eff += aS;
            }

            const double u_no_relax = rhs / aP_eff;
            u_star(i, j) = s.alpha_u * u_no_relax + (1.0 - s.alpha_u) * u_old(i, j);
            aP_u(i, j)   = aP_eff;
        }
    }
}

// ---------------------------------------------------------------------------
// Momentum predictor — v component
// Interior v faces: i in [0, Nx-1], j in [1, Ny-1]
// ---------------------------------------------------------------------------
static void v_momentum(const Mesh& m, const SolverSettings& s,
                        const Field& u_old, const Field& v_old,
                        const Field& p,
                        Field& v_star, Field& aP_v) {
    const double dx = m.dx, dy = m.dy;
    const double De = s.nu * dy / dx;
    const double Dn = s.nu * dx / dy;

    for (int j = 1; j < m.Ny; ++j) {
        for (int i = 0; i < m.Nx; ++i) {

            // u at east/west faces of the v-CV
            const double u_e = 0.5 * (u_old(i + 1, j - 1) + u_old(i + 1, j));
            const double u_w = 0.5 * (u_old(i,     j - 1) + u_old(i,     j));
            const double v_n = 0.5 * (v_old(i, j) + v_old(i, j + 1));
            const double v_s = 0.5 * (v_old(i, j - 1) + v_old(i, j));

            const double Fe = s.rho * u_e * dy;
            const double Fw = s.rho * u_w * dy;
            const double Fn = s.rho * v_n * dx;
            const double Fs = s.rho * v_s * dx;

            const double aE  = De + std::max(-Fe, 0.0);
            const double aW  = De + std::max( Fw, 0.0);
            const double aN  = Dn + std::max(-Fn, 0.0);
            const double aS  = Dn + std::max( Fs, 0.0);
            const double aP0 = s.rho * dx * dy / s.pseudo_dt;
            const double aP  = aE + aW + aN + aS + aP0;

            const double b = -(p(i, j) - p(i, j - 1)) * dx;

            double rhs    = aP0 * v_old(i, j) + b;
            rhs += aN * v_old(i, j + 1);  // top wall: v(i,Ny)=0, so index Ny is valid
            rhs += aS * v_old(i, j - 1);  // bottom wall: v(i,0)=0

            // East/west neighbours with inlet/outlet treatment
            double aP_eff = aP;
            if (i > 0) {
                rhs += aW * v_old(i - 1, j);
            } else {
                // Inlet: v=0 ghost gives ghost = -v_P → no-penetration
                aP_eff += aW;
            }
            if (i < m.Nx - 1) {
                rhs += aE * v_old(i + 1, j);
            } else {
                // Outlet: zero-gradient ghost = v_P → subtract from diagonal
                aP_eff -= aE;
            }

            const double v_no_relax = rhs / aP_eff;
            v_star(i, j) = s.alpha_v * v_no_relax + (1.0 - s.alpha_v) * v_old(i, j);
            aP_v(i, j)   = aP_eff;
        }
    }
}

// ---------------------------------------------------------------------------
// Pressure-correction Poisson equation — Conjugate Gradient solver
//
// System: A * pc = b, where b(i,j) = -div(u*) and A is the discrete
// Laplacian with coefficients aE = dy/(aP_u(i+1,j)*dx), etc.
//
// BCs:
//   Neumann at inlet (i=0 west) and walls (j=0/Ny south/north): omit face
//   Dirichlet at outlet: pc(Nx-1,j) = 0, treated as rhs contribution
//
// A is SPD (Dirichlet on east side pins the solution), so CG is exact.
// ---------------------------------------------------------------------------
static void pressure_correction(const Mesh& m, const SolverSettings& s,
                                 const Field& u_star, const Field& v_star,
                                 const Field& aP_u,   const Field& aP_v,
                                 Field& pc) {
    const double dx = m.dx, dy = m.dy;
    const int    Nx = m.Nx, Ny = m.Ny;
    const int    nu = (Nx - 1) * Ny;

    auto idx = [&](int i, int j) { return j * (Nx - 1) + i; };

    // RHS: b = -div(u*)
    std::vector<double> b(nu);
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx - 1; ++i)
            b[idx(i, j)] = -(u_star(i + 1, j) - u_star(i, j)) / dx
                           -(v_star(i, j + 1) - v_star(i, j)) / dy;

    // Matrix-vector product: y = A*x (A is SPD discrete Laplacian)
    auto matvec = [&](const std::vector<double>& x, std::vector<double>& y) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx - 1; ++i) {
                double diag = 0.0, off = 0.0;
                // East (always; Dirichlet pE=0 at outlet)
                { double a = dy / (aP_u(i + 1, j) * dx);
                  diag += a;
                  if (i < Nx - 2) off += a * x[idx(i + 1, j)]; }
                // West (Neumann at inlet: skip)
                if (i > 0) { double a = dy / (aP_u(i, j) * dx); diag += a; off += a * x[idx(i - 1, j)]; }
                // North (Neumann at wall: skip)
                if (j < Ny - 1) { double a = dx / (aP_v(i, j + 1) * dy); diag += a; off += a * x[idx(i, j + 1)]; }
                // South (Neumann at wall: skip)
                if (j > 0) { double a = dx / (aP_v(i, j) * dy); diag += a; off += a * x[idx(i, j - 1)]; }
                y[idx(i, j)] = diag * x[idx(i, j)] - off;
            }
        }
    };

    // CG: solve A*x = b starting from x=0
    const int    max_cg  = std::max(s.pressure_sweeps * 5, 500);
    const double tol_rel = 1.0e-6;

    std::vector<double> x(nu, 0.0);
    std::vector<double> r(b);          // r = b - A*0 = b
    std::vector<double> p(r);
    std::vector<double> Ap(nu);

    double rr = 0.0;
    for (int k = 0; k < nu; ++k) rr += r[k] * r[k];
    const double rr0 = rr;

    for (int it = 0; it < max_cg && rr > tol_rel * tol_rel * rr0; ++it) {
        matvec(p, Ap);
        double pAp = 0.0;
        for (int k = 0; k < nu; ++k) pAp += p[k] * Ap[k];
        const double alpha = rr / pAp;
        for (int k = 0; k < nu; ++k) { x[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }
        double rr_new = 0.0;
        for (int k = 0; k < nu; ++k) rr_new += r[k] * r[k];
        const double beta = rr_new / rr;
        rr = rr_new;
        for (int k = 0; k < nu; ++k) p[k] = r[k] + beta * p[k];
    }

    // Write solution back to pc field (outlet column stays 0)
    pc.fill(0.0);
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx - 1; ++i)
            pc(i, j) = x[idx(i, j)];
}

// ---------------------------------------------------------------------------
// Velocity and pressure correction
// ---------------------------------------------------------------------------
static void correct(const Mesh& m, const SolverSettings& s,
                    const Field& pc, const Field& aP_u, const Field& aP_v,
                    Field& u, Field& v, Field& p) {
    const double dx = m.dx, dy = m.dy;

    // u correction at interior faces
    for (int j = 0; j < m.Ny; ++j) {
        for (int i = 1; i < m.Nx; ++i) {
            const double pc_right = (i <= m.Nx - 1) ? pc(i, j) : 0.0;
            const double pc_left  = (i >= 1)         ? pc(i - 1, j) : pc(0, j);
            u(i, j) -= (dy / aP_u(i, j)) * (pc_right - pc_left);
        }
    }

    // v correction at interior faces
    for (int j = 1; j < m.Ny; ++j) {
        for (int i = 0; i < m.Nx; ++i) {
            const double pc_top = pc(i, j);
            const double pc_bot = pc(i, j - 1);
            v(i, j) -= (dx / aP_v(i, j)) * (pc_top - pc_bot);
        }
    }

    // Pressure update (under-relaxed)
    for (int j = 0; j < m.Ny; ++j) {
        for (int i = 0; i < m.Nx; ++i) {
            p(i, j) += s.alpha_p * pc(i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// RMS divergence of the corrected (u, v) on all pressure cells
// ---------------------------------------------------------------------------
// Reports two values: overall RMS and interior-only RMS (excluding outlet column)
static double divergence_residual(const Mesh& m, const Field& u, const Field& v,
                                  double* interior_div = nullptr) {
    double sum_all  = 0.0;
    double sum_int  = 0.0;
    int    cnt_int  = 0;
    for (int j = 0; j < m.Ny; ++j) {
        for (int i = 0; i < m.Nx; ++i) {
            const double div = (u(i + 1, j) - u(i, j)) / m.dx
                             + (v(i, j + 1) - v(i, j)) / m.dy;
            sum_all += div * div;
            if (i < m.Nx - 1) { sum_int += div * div; ++cnt_int; }
        }
    }
    if (interior_div) *interior_div = (cnt_int > 0) ? std::sqrt(sum_int / cnt_int) : 0.0;
    return std::sqrt(sum_all / (m.Nx * m.Ny));
}

// ---------------------------------------------------------------------------
// Helper: RMS change in a field (fluid cells only, here all cells)
// ---------------------------------------------------------------------------
static double rms_change(const Field& a, const Field& b) {
    double sum = 0.0;
    const int n = a.nx * a.ny;
    for (int k = 0; k < n; ++k) {
        const double d = a(k % a.nx, k / a.nx) - b(k % b.nx, k / b.nx);
        sum += d * d;
    }
    return std::sqrt(sum / n);
}

// ---------------------------------------------------------------------------
// Main SIMPLE loop
// ---------------------------------------------------------------------------
SolverReport run_simple(const Mesh& m, const SolverSettings& s, SolverState& st) {
    SolverReport report{};

    Field u_old(m.u_nx(), m.u_ny());
    Field v_old(m.v_nx(), m.v_ny());
    Field u_star(m.u_nx(), m.u_ny());
    Field v_star(m.v_nx(), m.v_ny());
    Field aP_u(m.u_nx(), m.u_ny(), 1.0);
    Field aP_v(m.v_nx(), m.v_ny(), 1.0);
    Field pc(m.p_nx(), m.p_ny());

    apply_velocity_bc(m, st.u, st.v, s.u_in);
    apply_pressure_bc(m, st.p);

    for (int iter = 1; iter <= s.max_outer; ++iter) {
        u_old = st.u;
        v_old = st.v;

        // --- predictor ---
        u_star = u_old;
        v_star = v_old;
        u_momentum(m, s, u_old, v_old, st.p, u_star, aP_u);
        v_momentum(m, s, u_old, v_old, st.p, v_star, aP_v);
        apply_velocity_bc(m, u_star, v_star, s.u_in);

        // --- pressure correction ---
        pressure_correction(m, s, u_star, v_star, aP_u, aP_v, pc);

        // --- corrector ---
        st.u = u_star;
        st.v = v_star;
        correct(m, s, pc, aP_u, aP_v, st.u, st.v, st.p);
        apply_velocity_bc(m, st.u, st.v, s.u_in);
        apply_pressure_bc(m, st.p);

        // Divergence of the corrected velocity (convergence criterion)
        double div_int = 0.0;
        report.div_residual = divergence_residual(m, st.u, st.v, &div_int);

        if (iter % s.report_interval == 0 || iter == 1) {
            std::cout << "  [div_interior=" << div_int << "]\n";
        }

        report.iterations = iter;
        report.u_change   = rms_change(st.u, u_old);
        report.v_change   = rms_change(st.v, v_old);

        if (iter % s.report_interval == 0 || iter == 1) {
            std::cout << "iter " << iter
                      << "  div=" << report.div_residual
                      << "  du=" << report.u_change
                      << "  dv=" << report.v_change << '\n';
        }

        if (report.u_change < s.tol_u) {
            std::cout << "Converged at iter " << iter
                      << "  du=" << report.u_change
                      << "  div=" << report.div_residual << '\n';
            break;
        }
    }

    return report;
}
