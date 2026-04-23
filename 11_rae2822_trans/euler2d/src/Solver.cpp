#include "Solver.hpp"
#include "BC.hpp"
#include "Flux.hpp"
#include "IO.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

void Solver::set_cons(int i, int j, double r, double u, double v, double p) {
    U(0, i, j) = r;
    U(1, i, j) = r * u;
    U(2, i, j) = r * v;
    U(3, i, j) = p / (gamma - 1.0) + 0.5 * r * (u*u + v*v);
}

double Solver::press(int i, int j) const {
    double r  = U(0, i, j);
    double ru = U(1, i, j);
    double rv = U(2, i, j);
    double rE = U(3, i, j);
    return (gamma - 1.0) * (rE - 0.5 * (ru*ru + rv*rv) / r);
}

double Solver::sound(int i, int j) const {
    return std::sqrt(gamma * press(i, j) / U(0, i, j));
}

void Solver::init(const Mesh& m,
                  double mach_, double aoa_deg_, double gamma_,
                  double cfl_, int max_iter_, double residual_drop_,
                  int scheme_order_, int output_interval_) {
    pmesh           = &m;
    gamma           = gamma_;
    mach            = mach_;
    aoa_deg         = aoa_deg_;
    cfl             = cfl_;
    max_iter        = max_iter_;
    residual_drop   = residual_drop_;
    scheme_order    = scheme_order_;
    output_interval = output_interval_;

    nci = m.nc_i();
    ncj = m.nc_j();

    rho_inf = 1.0;
    p_inf   = 1.0 / gamma;
    double aoa_rad = aoa_deg * std::acos(-1.0) / 180.0;
    u_inf   = mach * std::cos(aoa_rad);
    v_inf   = mach * std::sin(aoa_rad);
    E_inf   = p_inf / (gamma - 1.0) + 0.5 * (u_inf*u_inf + v_inf*v_inf);

    data.assign(4 * (nci + 2) * (ncj + 2), 0.0);
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            set_cons(i, j, rho_inf, u_inf, v_inf, p_inf);
    apply_bc();
}

void Solver::apply_bc() {
    bc_wall(*this);
    bc_farfield(*this);
    bc_wake_cut(*this);
    bc_i_exit(*this);
}

void Solver::print_bc_diagnostics() const {
    int iw = nci / 2;
    std::cout << "\n=== BC diagnostics ===\n";
    std::cout << "Wall ghost (i=" << iw << ", j=-1):  rho=" << rho(iw,-1)
              << " u=" << vel_u(iw,-1) << " v=" << vel_v(iw,-1)
              << " p=" << press(iw,-1) << "\n";
    std::cout << "Wall real  (i=" << iw << ", j=0):   rho=" << rho(iw,0)
              << " u=" << vel_u(iw,0)  << " v=" << vel_v(iw,0)
              << " p=" << press(iw,0)  << "\n";
    std::cout << "Farfield ghost (i=" << iw << ", j=" << ncj << "): rho=" << rho(iw,ncj)
              << " u=" << vel_u(iw,ncj) << " v=" << vel_v(iw,ncj)
              << " p=" << press(iw,ncj) << "\n";
    std::cout << "Freestream: rho=" << rho_inf
              << " u=" << u_inf << " v=" << v_inf << " p=" << p_inf << "\n";
    int i_wc = 30, i_mirror = nci - 1 - i_wc;
    std::cout << "Wake-cut ghost (i=" << i_wc << ",j=-1) rho=" << rho(i_wc,-1)
              << "  mirror (i=" << i_mirror << ",j=0) rho=" << rho(i_mirror,0) << "\n";
}

// ---- compute_rhs ------------------------------------------------------------
// Fills R[var * N + j*nci + i] (N=nci*ncj) and dt_cell[j*nci+i].
// Sign convention: dU/dt = R  →  U += dt * R.
// For each face with outward normal n (from left cell L to right cell R):
//   R[L] -= F_roe * len / area[L]
//   R[R] += F_roe * len / area[R]

void Solver::compute_rhs(std::vector<double>& R, std::vector<double>& dt_cell) const {
    const Mesh& m = *pmesh;
    int N  = nci * ncj;
    int ni = m.ni;   // node count in i = nci+1

    std::fill(R.begin(), R.end(), 0.0);

    // Spectral-radius accumulator per cell for local time step
    std::vector<double> srad(N, 0.0);

    // Helper: get primitives at cell (ci, cj) — supports ghost
    auto prim = [&](int ci, int cj, double& r, double& u, double& v, double& p) {
        r = rho  (ci, cj);
        u = vel_u(ci, cj);
        v = vel_v(ci, cj);
        p = press(ci, cj);
    };

    // ---- i-direction faces: face i is between cell (i-1,j) [L] and (i,j) [R]
    //      i loops 0..nci  (ni = nci+1 faces total)
    //      ifnx/ifny indexed [j*ni + i]
    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i <= nci; ++i) {
            int iL = i - 1, iR = i;
            double nx  = m.ifnx[j * ni + i];
            double ny  = m.ifny[j * ni + i];
            double len = m.iflen[j * ni + i];
            if (len < 1e-14) continue;  // skip degenerate faces (wake-exit singularity)

            double rL, uL, vL, pL, rR, uR, vR, pR;
            prim(iL, j, rL, uL, vL, pL);
            prim(iR, j, rR, uR, vR, pR);

            double F[4];
            roe_flux(rL, uL, vL, pL, rR, uR, vR, pR, nx, ny, gamma, F);

            // Spectral radius for local dt
            double cL = std::sqrt(gamma * pL / rL);
            double cR = std::sqrt(gamma * pR / rR);
            double sr = (std::abs(uL*nx + vL*ny) + cL +
                         std::abs(uR*nx + vR*ny) + cR) * 0.5 * len;

            if (iL >= 0 && iL < nci) {
                int idx = j * nci + iL;
                double invA = 1.0 / m.area[idx];
                for (int v = 0; v < 4; ++v) R[v*N + idx] -= F[v] * len * invA;
                srad[idx] += sr;
            }
            if (iR >= 0 && iR < nci) {
                int idx = j * nci + iR;
                double invA = 1.0 / m.area[idx];
                for (int v = 0; v < 4; ++v) R[v*N + idx] += F[v] * len * invA;
                srad[idx] += sr;
            }
        }
    }

    // ---- j-direction faces: face j is between cell (i,j-1) [L] and (i,j) [R]
    //      j loops 0..ncj  (nj = ncj+1 faces total)
    //      jfnx/jfny indexed [j*nci + i]
    for (int j = 0; j <= ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
            int jL = j - 1, jR = j;
            double nx  = m.jfnx[j * nci + i];
            double ny  = m.jfny[j * nci + i];
            double len = m.jflen[j * nci + i];
            if (len < 1e-14) continue;  // skip degenerate faces

            double rL, uL, vL, pL, rR, uR, vR, pR;
            prim(i, jL, rL, uL, vL, pL);
            prim(i, jR, rR, uR, vR, pR);

            double F[4];
            roe_flux(rL, uL, vL, pL, rR, uR, vR, pR, nx, ny, gamma, F);

            double cL = std::sqrt(gamma * pL / rL);
            double cR = std::sqrt(gamma * pR / rR);
            double sr = (std::abs(uL*nx + vL*ny) + cL +
                         std::abs(uR*nx + vR*ny) + cR) * 0.5 * len;

            if (jL >= 0 && jL < ncj) {
                int idx = jL * nci + i;
                double invA = 1.0 / m.area[idx];
                for (int v = 0; v < 4; ++v) R[v*N + idx] -= F[v] * len * invA;
                srad[idx] += sr;
            }
            if (jR >= 0 && jR < ncj) {
                int idx = jR * nci + i;
                double invA = 1.0 / m.area[idx];
                for (int v = 0; v < 4; ++v) R[v*N + idx] += F[v] * len * invA;
                srad[idx] += sr;
            }
        }
    }

    // Local time steps
    for (int k = 0; k < N; ++k)
        dt_cell[k] = (srad[k] > 1e-30) ? cfl * m.area[k] / srad[k] : 0.0;

}

// ---- SSP-RK3 run loop -------------------------------------------------------

void Solver::run() {
    int N = nci * ncj;
    std::vector<double> R(4 * N, 0.0);
    std::vector<double> dt_cell(N, 0.0);

    // Backup of real-cell state at start of each RK step
    // Stored as Un[var*N + j*nci + i]
    std::vector<double> Un(4 * N, 0.0);

    auto pack_real = [&](std::vector<double>& buf) {
        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int v = 0; v < 4; ++v)
                    buf[v*N + j*nci + i] = U(v, i, j);
    };
    std::ofstream conv("output/convergence.dat");
    conv << "# iter  res_rho  res_norm\n";

    double res0 = -1.0;

    for (int iter = 1; iter <= max_iter; ++iter) {

        // Save Un for RK combination
        pack_real(Un);

        // ---- RK stage 1: U1 = Un + dt*R(Un) ----
        apply_bc();
        compute_rhs(R, dt_cell);

        // Convergence metric (density residual from stage 1)
        double res_rho = 0.0;
        for (int k = 0; k < N; ++k) res_rho += R[k] * R[k];
        res_rho = std::sqrt(res_rho / N);
        if (iter == 1) res0 = res_rho;
        double res_norm = (res0 > 0.0) ? res_rho / res0 : 1.0;

        conv << iter << " " << res_rho << " " << res_norm << "\n";
        if (iter % output_interval == 0 || iter <= 2)
            std::cout << "iter " << iter
                      << "  res=" << res_rho
                      << "  res/res0=" << res_norm << "\n";

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int v = 0; v < 4; ++v)
                    U(v, i, j) = Un[v*N + j*nci+i] + dt_cell[j*nci+i] * R[v*N + j*nci+i];

        // ---- RK stage 2: U2 = 0.75*Un + 0.25*(U1 + dt*R(U1)) ----
        apply_bc();
        compute_rhs(R, dt_cell);

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int v = 0; v < 4; ++v) {
                    double u1 = U(v, i, j);
                    double un = Un[v*N + j*nci+i];
                    U(v, i, j) = 0.75*un + 0.25*(u1 + dt_cell[j*nci+i]*R[v*N + j*nci+i]);
                }

        // ---- RK stage 3: Un+1 = (1/3)*Un + (2/3)*(U2 + dt*R(U2)) ----
        apply_bc();
        compute_rhs(R, dt_cell);

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int v = 0; v < 4; ++v) {
                    double u2 = U(v, i, j);
                    double un = Un[v*N + j*nci+i];
                    U(v, i, j) = (1.0/3.0)*un + (2.0/3.0)*(u2 + dt_cell[j*nci+i]*R[v*N + j*nci+i]);
                }

        // Write output files periodically and at convergence
        bool last = (res_norm < residual_drop) || (iter == max_iter);
        if (iter % output_interval == 0 || last) {
            write_field(*this, "output/field.vts");
            write_cp   (*this, "output/cp.dat");
        }

        if (res_norm < residual_drop) {
            std::cout << "Converged at iter " << iter
                      << "  res/res0=" << res_norm << "\n";
            break;
        }
    }
}
