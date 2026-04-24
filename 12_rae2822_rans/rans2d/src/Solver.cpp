#include "Solver.hpp"
#include "BC.hpp"
#include "Flux.hpp"
#include "IO.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

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
                  double cfl_, double cfl_impl_, double omega_,
                  int max_iter_, double residual_drop_,
                  int scheme_order_, int warmup_iters_,
                  int output_interval_) {
    pmesh         = &m;
    gamma         = gamma_;
    mach          = mach_;
    aoa_deg       = aoa_deg_;
    cfl           = cfl_;
    cfl_impl      = cfl_impl_;
    omega         = omega_;
    max_iter      = max_iter_;
    residual_drop = residual_drop_;
    scheme_order  = scheme_order_;
    warmup_iters  = warmup_iters_;
    output_interval = output_interval_;

    nci = m.nc_i();
    ncj = m.nc_j();

    rho_inf = 1.0;
    p_inf   = 1.0 / gamma;
    double aoa_rad = aoa_deg * std::acos(-1.0) / 180.0;
    u_inf = mach * std::cos(aoa_rad);
    v_inf = mach * std::sin(aoa_rad);
    E_inf = p_inf / (gamma - 1.0) + 0.5 * (u_inf*u_inf + v_inf*v_inf);

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
}

// Cells whose area falls below this threshold are "bridge cells" — topological
// artefacts of the C-grid where two wall nodes share the same physical point.
// These cells are frozen (R=0, dU=0) to prevent residual blow-up from invA→∞.
// TE singularity: column i=191 has nodes coinciding at (1,0) for j=0.
// Cells in this column have areas that grow from ~6e-14 upward.
// Standard first-layer airfoil cells have area ~3e-10 c².
// Freeze all cells below a safe margin below the first real cell.
static constexpr double AREA_MIN = 1e-9;

// ---- compute_rhs ------------------------------------------------------------
// Fills R[var*N + j*nci + i] (N=nci*ncj) and dt_cell[j*nci+i].
// Sign: dU/dt = R.  Degenerate bridge cells (area < AREA_MIN) get R=0, dt=0.

void Solver::compute_rhs(std::vector<double>& R,
                          std::vector<double>& dt_cell) const {
    const Mesh& m = *pmesh;
    int N  = nci * ncj;
    int ni = m.ni;

    std::fill(R.begin(), R.end(), 0.0);
    std::vector<double> srad_cell(N, 0.0);

    auto prim = [&](int ci, int cj, double& r, double& u, double& v, double& p) {
        r = rho  (ci, cj);
        u = vel_u(ci, cj);
        v = vel_v(ci, cj);
        p = press(ci, cj);
    };

    auto muscl = [](double bm1, double b0, double b1, double b2,
                    double mf, double& qL, double& qR) {
        auto phi = [](double r) -> double {
            return r > 0.0 ? (r*r + r) / (r*r + 1.0) : 0.0;  // van Albada
        };
        double dl = b0 - bm1, dc = b1 - b0, dr = b2 - b1;
        double rL = (std::abs(dl) > 1e-14) ? dc / dl : 0.0;
        double rR = (std::abs(dr) > 1e-14) ? dc / dr : 0.0;
        qL = b0 + 0.5 * mf * phi(rL) * dl;
        qR = b1 - 0.5 * mf * phi(rR) * dr;
    };

    // Pressure-based shock sensor: 0 at shock, 1 in smooth flow.
    auto shock_alpha = [](double pA, double p0, double p1, double pB) -> double {
        const double eps = 1e-10;
        double s = std::max({std::abs(p0-pA)/(p0+pA+eps),
                             std::abs(p1-p0)/(p1+p0+eps),
                             std::abs(pB-p1)/(pB+p1+eps)});
        return std::max(0.0, 1.0 - 3.0 * s);
    };

    // ---- i-direction faces ----
    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i <= nci; ++i) {
            int iL = i - 1, iR = i;
            double nx  = m.ifnx[j*ni + i];
            double ny  = m.ifny[j*ni + i];
            double len = m.iflen[j*ni + i];
            if (len < 1e-14) continue;

            double rho0, u0, v0, p0, rho1, u1, v1, p1;
            prim(iL, j, rho0, u0, v0, p0);
            prim(iR, j, rho1, u1, v1, p1);

            double rL, uL, vL, pL, rR, uR, vR, pR;
            // Disable MUSCL if any REAL stencil cell is a frozen bridge cell.
            // Ghost cells (index outside [0,nci)) don't have areas — skip them.
            bool ok_muscl_i = (scheme_order >= 2 && i >= 1 && i <= nci - 1)
                && m.area[j*nci+iL] >= AREA_MIN
                && m.area[j*nci+iR] >= AREA_MIN
                && (iL-1 < 0    || m.area[j*nci+(iL-1)] >= AREA_MIN)
                && (iR+1 >= nci || m.area[j*nci+(iR+1)] >= AREA_MIN);
            if (ok_muscl_i) {
                double rhoA, uA, vA, pA, rhoB, uB, vB, pB;
                prim(iL-1, j, rhoA, uA, vA, pA);
                prim(iR+1, j, rhoB, uB, vB, pB);
                double mf = muscl_factor * shock_alpha(pA, p0, p1, pB);
                muscl(rhoA, rho0, rho1, rhoB, mf, rL, rR);
                muscl(uA,   u0,   u1,   uB,   mf, uL, uR);
                muscl(vA,   v0,   v1,   vB,   mf, vL, vR);
                muscl(pA,   p0,   p1,   pB,   mf, pL, pR);
            } else {
                rL=rho0; uL=u0; vL=v0; pL=p0;
                rR=rho1; uR=u1; vR=v1; pR=p1;
            }
            rL=std::max(rL,1e-10); rR=std::max(rR,1e-10);
            pL=std::max(pL,1e-10); pR=std::max(pR,1e-10);

            double F[4];
            roe_flux(rL, uL, vL, pL, rR, uR, vR, pR, nx, ny, gamma, F);

            double cL = std::sqrt(gamma*pL/rL);
            double cR = std::sqrt(gamma*pR/rR);
            double sr = (std::abs(uL*nx+vL*ny)+cL + std::abs(uR*nx+vR*ny)+cR) * 0.5 * len;

            if (iL >= 0 && iL < nci) {
                int idx = j*nci + iL;
                if (m.area[idx] >= AREA_MIN) {
                    double invA = 1.0 / m.area[idx];
                    for (int v = 0; v < 4; ++v) R[v*N + idx] -= F[v] * len * invA;
                    srad_cell[idx] += sr;
                }
            }
            if (iR >= 0 && iR < nci) {
                int idx = j*nci + iR;
                if (m.area[idx] >= AREA_MIN) {
                    double invA = 1.0 / m.area[idx];
                    for (int v = 0; v < 4; ++v) R[v*N + idx] += F[v] * len * invA;
                    srad_cell[idx] += sr;
                }
            }
        }
    }

    // ---- j-direction faces ----
    for (int j = 0; j <= ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
            int jL = j-1, jR = j;
            double nx  = m.jfnx[j*nci + i];
            double ny  = m.jfny[j*nci + i];
            double len = m.jflen[j*nci + i];
            if (len < 1e-14) continue;

            double rho0, u0, v0, p0, rho1, u1, v1, p1;
            prim(i, jL, rho0, u0, v0, p0);
            prim(i, jR, rho1, u1, v1, p1);

            double rL, uL, vL, pL, rR, uR, vR, pR;
            // Disable j-direction MUSCL: wake-cut ghost creates artificial
            // gradients at j=1 faces that cause a growing limit cycle.
            // i-direction MUSCL alone is sufficient to capture the shock.
            bool ok_muscl_j = false;
            if (ok_muscl_j) {
                double rhoA, uA, vA, pA, rhoB, uB, vB, pB;
                prim(i, jL-1, rhoA, uA, vA, pA);
                prim(i, jR+1, rhoB, uB, vB, pB);
                double mf = muscl_factor * shock_alpha(pA, p0, p1, pB);
                muscl(rhoA, rho0, rho1, rhoB, mf, rL, rR);
                muscl(uA,   u0,   u1,   uB,   mf, uL, uR);
                muscl(vA,   v0,   v1,   vB,   mf, vL, vR);
                muscl(pA,   p0,   p1,   pB,   mf, pL, pR);
            } else {
                rL=rho0; uL=u0; vL=v0; pL=p0;
                rR=rho1; uR=u1; vR=v1; pR=p1;
            }
            rL=std::max(rL,1e-10); rR=std::max(rR,1e-10);
            pL=std::max(pL,1e-10); pR=std::max(pR,1e-10);

            double F[4];
            roe_flux(rL, uL, vL, pL, rR, uR, vR, pR, nx, ny, gamma, F);

            double cL = std::sqrt(gamma*pL/rL);
            double cR = std::sqrt(gamma*pR/rR);
            double sr = (std::abs(uL*nx+vL*ny)+cL + std::abs(uR*nx+vR*ny)+cR) * 0.5 * len;

            if (jL >= 0 && jL < ncj) {
                int idx = jL*nci + i;
                if (m.area[idx] >= AREA_MIN) {
                    double invA = 1.0 / m.area[idx];
                    for (int v = 0; v < 4; ++v) R[v*N + idx] -= F[v] * len * invA;
                    srad_cell[idx] += sr;
                }
            }
            if (jR >= 0 && jR < ncj) {
                int idx = jR*nci + i;
                if (m.area[idx] >= AREA_MIN) {
                    double invA = 1.0 / m.area[idx];
                    for (int v = 0; v < 4; ++v) R[v*N + idx] += F[v] * len * invA;
                    srad_cell[idx] += sr;
                }
            }
        }
    }

    // Local time steps from spectral radius sum
    for (int k = 0; k < N; ++k)
        dt_cell[k] = (srad_cell[k] > 1e-30) ? cfl * m.area[k] / srad_cell[k] : 0.0;
}

// ---- LU-SGS -----------------------------------------------------------------
//
// Solves (D+L+U) ΔU = area*R using the factored forward/backward sweep.
//
// Spectral radii at each face use the current cell's state ("frozen coeff"):
//   sr = (|u·n| + c) * len   (units: area/time)
//
// Diagonal: D = srad * (0.5 + 1/CFL)
//   = srad/CFL (explicit CFL term) + 0.5*srad (implicit spectral radius)
//
// Forward:   D * dU*_k = area*R_k + 0.5 * Σ_{l<k} sr_kl * dU*_l
// Backward:  dU_k = dU*_k + (0.5/D) * Σ_{l>k} sr_kl * dU_l

void Solver::lusgs(const std::vector<double>& R, std::vector<double>& dU) const {
    const Mesh& m = *pmesh;
    int N  = nci * ncj;
    int ni = m.ni;

    std::vector<double> dU_fwd(4*N, 0.0);
    std::vector<double> diag(N);
    std::vector<double> srR_store(N), srT_store(N);

    // Forward sweep: j = 0..ncj-1, i = 0..nci-1
    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
            int k = j*nci + i;
            double u  = vel_u(i,j), v  = vel_v(i,j), c = sound(i,j);
            double area = m.area[k];

            double nxL = m.ifnx[j*ni+i],     nyL = m.ifny[j*ni+i];
            double nxR = m.ifnx[j*ni+i+1],   nyR = m.ifny[j*ni+i+1];
            double nxB = m.jfnx[j*nci+i],    nyB = m.jfny[j*nci+i];
            double nxT = m.jfnx[(j+1)*nci+i],nyT = m.jfny[(j+1)*nci+i];

            double srL = (std::abs(u*nxL+v*nyL) + c) * m.iflen[j*ni+i];
            double srR = (std::abs(u*nxR+v*nyR) + c) * m.iflen[j*ni+i+1];
            double srB = (std::abs(u*nxB+v*nyB) + c) * m.jflen[j*nci+i];
            double srT = (std::abs(u*nxT+v*nyT) + c) * m.jflen[(j+1)*nci+i];

            double srad = srL + srR + srB + srT;
            double D = (srad > 1e-30) ? srad * (0.5 + 1.0/cfl) : 1.0;

            // Freeze degenerate bridge cells (near-zero area, C-grid TE artefact)
            if (area < AREA_MIN) {
                diag[k] = 1e30;  // coeff→0 in backward sweep
                srR_store[k] = 0.0;
                srT_store[k] = 0.0;
                for (int vv = 0; vv < 4; ++vv) dU_fwd[vv*N + k] = 0.0;
                continue;
            }

            diag[k]      = D;
            srR_store[k] = srR;
            srT_store[k] = srT;

            double rhs[4];
            for (int vv = 0; vv < 4; ++vv)
                rhs[vv] = R[vv*N + k] * area;

            // Off-diagonal from lower neighbors, using NEIGHBOR's spectral radius
            // at the shared face (Yoon-Jameson formulation for better stability on
            // high-AR grids — avoids frozen-coefficient over-amplification).
            if (i > 0) {
                int kl = j*nci + (i-1);
                double ul = vel_u(i-1,j), vl = vel_v(i-1,j), cl = sound(i-1,j);
                double srL_neigh = (std::abs(ul*nxL+vl*nyL) + cl) * m.iflen[j*ni+i];
                for (int vv = 0; vv < 4; ++vv)
                    rhs[vv] += 0.5 * srL_neigh * dU_fwd[vv*N + kl];
            }
            if (j > 0) {
                int kl = (j-1)*nci + i;
                double ul = vel_u(i,j-1), vl = vel_v(i,j-1), cl = sound(i,j-1);
                double nxB_neigh = m.jfnx[j*nci+i], nyB_neigh = m.jfny[j*nci+i];
                double srB_neigh = (std::abs(ul*nxB_neigh+vl*nyB_neigh) + cl) * m.jflen[j*nci+i];
                for (int vv = 0; vv < 4; ++vv)
                    rhs[vv] += 0.5 * srB_neigh * dU_fwd[vv*N + kl];
            }

            for (int vv = 0; vv < 4; ++vv)
                dU_fwd[vv*N + k] = rhs[vv] / D;
        }
    }

    // Backward sweep: j = ncj-1..0, i = nci-1..0
    // Use upper neighbor's spectral radius at the shared face.
    for (int j = ncj-1; j >= 0; --j) {
        for (int i = nci-1; i >= 0; --i) {
            int k = j*nci + i;
            double D = diag[k];

            if (m.area[k] < AREA_MIN) {
                for (int vv = 0; vv < 4; ++vv) dU[vv*N + k] = 0.0;
                continue;
            }

            for (int vv = 0; vv < 4; ++vv)
                dU[vv*N + k] = dU_fwd[vv*N + k];

            if (i < nci-1) {
                int ku = j*nci + (i+1);
                double uu = vel_u(i+1,j), vu = vel_v(i+1,j), cu = sound(i+1,j);
                double nxR_u = m.ifnx[j*ni+i+1], nyR_u = m.ifny[j*ni+i+1];
                double srR_u = (std::abs(uu*nxR_u+vu*nyR_u) + cu) * m.iflen[j*ni+i+1];
                double coeff = 0.5 * srR_u / D;
                for (int vv = 0; vv < 4; ++vv)
                    dU[vv*N + k] += coeff * dU[vv*N + ku];
            }
            if (j < ncj-1) {
                int ku = (j+1)*nci + i;
                double uu = vel_u(i,j+1), vu = vel_v(i,j+1), cu = sound(i,j+1);
                double nxT_u = m.jfnx[(j+1)*nci+i], nyT_u = m.jfny[(j+1)*nci+i];
                double srT_u = (std::abs(uu*nxT_u+vu*nyT_u) + cu) * m.jflen[(j+1)*nci+i];
                double coeff = 0.5 * srT_u / D;
                for (int vv = 0; vv < 4; ++vv)
                    dU[vv*N + k] += coeff * dU[vv*N + ku];
            }
        }
    }
}

// ---- run loop (SSP-RK3 + local time stepping) --------------------------------
// Same formulation as euler2d/src/Solver.cpp but with bridge-cell aware R.

void Solver::run() {
    int N = nci * ncj;
    std::vector<double> R(4*N, 0.0);
    std::vector<double> dt_cell(N, 0.0);
    std::vector<double> Un(4*N, 0.0);

    int final_order = scheme_order;
    if (warmup_iters > 0 && scheme_order >= 2) {
        scheme_order = 1;
        std::cout << "  [warmup] order=1 for " << warmup_iters << " iters\n";
    }

    auto pack_real = [&](std::vector<double>& buf) {
        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int vv = 0; vv < 4; ++vv)
                    buf[vv*N + j*nci+i] = U(vv, i, j);
    };

    std::ofstream conv("output/convergence.dat");
    conv << "# iter  res_rho  res_norm\n";
    double res0 = -1.0;

    for (int iter = 1; iter <= max_iter; ++iter) {

        if (scheme_order < final_order && iter > warmup_iters) {
            scheme_order = final_order;
            cfl          = cfl_impl;
            res0         = -1.0;
            std::cout << "  [order switch] order=" << scheme_order
                      << "  CFL=" << cfl << "  at iter " << iter << "\n";
            conv << "# order switch iter=" << iter << "\n";
        }
        muscl_factor = (scheme_order >= 2 && iter > warmup_iters) ? 1.0 : 0.0;

        pack_real(Un);

        // ---- RK stage 1 ----
        apply_bc();
        compute_rhs(R, dt_cell);

        double res_rho = 0.0;
        for (int k = 0; k < N; ++k) res_rho += R[k] * R[k];
        res_rho = std::sqrt(res_rho / N);
        if (res0 < 0.0) res0 = res_rho;
        double res_norm = (res0 > 0.0) ? res_rho / res0 : 1.0;

        conv << iter << " " << res_rho << " " << res_norm << "\n";
        if (iter % output_interval == 0 || iter <= 2 || iter == warmup_iters + 1)
            std::cout << "iter " << std::setw(6) << iter
                      << "  res=" << res_rho
                      << "  res/res0=" << res_norm << "\n";

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int vv = 0; vv < 4; ++vv)
                    U(vv,i,j) = Un[vv*N + j*nci+i] + dt_cell[j*nci+i]*R[vv*N + j*nci+i];

        // ---- RK stage 2 ----
        apply_bc();
        compute_rhs(R, dt_cell);

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int vv = 0; vv < 4; ++vv) {
                    double u1 = U(vv,i,j);
                    double un = Un[vv*N + j*nci+i];
                    U(vv,i,j) = 0.75*un + 0.25*(u1 + dt_cell[j*nci+i]*R[vv*N + j*nci+i]);
                }

        // ---- RK stage 3 ----
        apply_bc();
        compute_rhs(R, dt_cell);

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int vv = 0; vv < 4; ++vv) {
                    double u2 = U(vv,i,j);
                    double un = Un[vv*N + j*nci+i];
                    U(vv,i,j) = (1.0/3.0)*un + (2.0/3.0)*(u2 + dt_cell[j*nci+i]*R[vv*N + j*nci+i]);
                }

        bool last = (res_norm < residual_drop) || (iter == max_iter);
        if (iter % output_interval == 0 || last) {
            write_field  (*this, "output/field.vts");
            write_surface(*this, "output/surface.vtp");
            write_cp     (*this, "output/cp.dat");
        }

        if (res_norm < residual_drop) {
            std::cout << "Converged at iter " << iter
                      << "  res/res0=" << res_norm << "\n";
            break;
        }
    }
}
