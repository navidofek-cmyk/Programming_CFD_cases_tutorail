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

// ---- Spalart-Allmaras helpers -----------------------------------------------

double Solver::sa_fv1(double chi) const {
    double c3 = SA_cv1 * SA_cv1 * SA_cv1;
    double x3 = chi * chi * chi;
    return x3 / (x3 + c3);
}

double Solver::sa_fw(double r) const {
    double g  = r + SA_cw2 * (std::pow(r, 6.0) - r);
    double c6 = std::pow(SA_cw3, 6.0);
    return g * std::pow((1.0 + c6) / (std::pow(g, 6.0) + c6), 1.0 / 6.0);
}

double Solver::mu_t(int i, int j) const {
    if (reynolds <= 0.0 || sa_data.empty()) return 0.0;
    double chi  = SA(i,j) * reynolds / rho(i,j);  // ν̃ / ν_mol = ρν̃*Re/ρ
    return SA(i,j) * sa_fv1(chi);                  // ρν̃ * fv1  (= ρ*ν_t, non-dim)
}

double Solver::mu_eff(int i, int j) const {
    // μ_eff = μ + μ_t = 1/Re + ρν̃*fv1
    double mu_lam = (reynolds > 0.0) ? 1.0 / reynolds : 0.0;
    return mu_lam + mu_t(i, j);
}

double Solver::kappa_eff(int i, int j) const {
    // κ_eff: lam + turb Prandtl contributions to heat conductivity
    // κ = μ * γ / ((γ-1) * Pr)
    double mu_lam = (reynolds > 0.0) ? 1.0 / reynolds : 0.0;
    return mu_lam * gamma / ((gamma - 1.0) * prandtl)
         + mu_t(i, j) * gamma / ((gamma - 1.0) * prandtl_t);
}

void Solver::init(const Mesh& m,
                  double mach_, double aoa_deg_, double gamma_,
                  double reynolds_, double prandtl_, double prandtl_t_,
                  double cfl_, double cfl_impl_, double omega_,
                  int max_iter_, double residual_drop_,
                  int scheme_order_, int warmup_iters_,
                  int output_interval_) {
    pmesh         = &m;
    gamma         = gamma_;
    mach          = mach_;
    aoa_deg       = aoa_deg_;
    reynolds      = reynolds_;
    prandtl       = prandtl_;
    prandtl_t     = prandtl_t_;
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

    // SA initialization: start with ν̃ = 0 (cold start).
    // The freestream BC maintains ñ_inf=3/Re at the farfield.
    // Starting from 0 avoids the large initial destruction term that arises
    // when S̃ < 0 in uniform flow (no wall shear to sustain turbulence).
    double nu_inf = (reynolds > 0.0) ? 3.0 / reynolds : 0.0;
    // SA initialization: ν̃ = 3 * ν_mol = 3/Re throughout domain.
    // Semi-implicit destruction prevents stiff blow-up during the initial
    // transient where S̃ < 0 (no shear in freestream).
    sa_data.assign((nci + 2) * (ncj + 2), 0.0);
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            SA(i, j) = rho_inf * nu_inf;

    apply_bc();
}

void Solver::apply_bc() {
    if (reynolds > 0.0) bc_wall_noslip(*this);
    else                bc_wall(*this);
    bc_farfield(*this);
    bc_wake_cut(*this);
    bc_i_exit(*this);

    if (sa_data.empty()) return;  // SA not initialised yet

    // ---- SA ghost cells ----
    double nu_inf  = (reynolds > 0.0) ? 3.0 / reynolds : 0.0;
    double sa_inf  = rho_inf * nu_inf;   // ρν̃ at farfield

    // Wall j=-1: ν̃=0 → ρν̃_ghost = -ρν̃_real (mirror)
    for (int i = pmesh->i_TEl; i < pmesh->i_TEu; ++i)
        SA(i, -1) = -SA(i, 0);

    // Farfield j=ncj: ρν̃ = sa_inf
    for (int i = 0; i < nci; ++i)
        SA(i, ncj) = sa_inf;

    // Wake cut: mirror (same as NS variables)
    for (int i = 0; i < pmesh->i_TEl; ++i)
        SA(i, -1) = SA(nci - 1 - i, 0);
    for (int i = pmesh->i_TEu; i < nci; ++i)
        SA(i, -1) = SA(nci - 1 - i, 0);

    // i-exit: extrapolation
    for (int j = 0; j < ncj; ++j) {
        SA(-1,  j) = SA(0,      j);
        SA(nci, j) = SA(nci-1,  j);
    }
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

    // ---- Viscous fluxes (Stage 3+) — Thin-Layer NS ----------------------------
    // Only when reynolds > 0.  Thin-layer approximation: only j-direction (wall-
    // normal) viscous fluxes.  Gradients computed by direct cell-center difference
    // rather than Green-Gauss, which is unstable on highly-stretched RANS meshes.
    //
    // Dominant physics: ∂τ/∂n drives wall friction; i-direction viscous terms are
    // O(1/AR) smaller and omitted (thin-layer assumption, valid for Re >> 1).
    //
    // Non-dimensional: μ = 1/Re, κ = μ·γ/((γ-1)·Pr)
    // Signs: R[L] += Fv·len/area_L, R[R] -= Fv·len/area_R

    if (reynolds <= 0.0) return;

    // Helper: build Fv from normal gradient, with per-face mu_f and kappa_f
    auto viscous_face = [&](int kL, int kR,
                            double nx, double ny, double len2,
                            double dudn, double dvdn, double dTdn,
                            double uf, double vf,
                            double mu_f, double kappa_f) {
        double dudx = dudn*nx, dudy = dudn*ny;
        double dvdx = dvdn*nx, dvdy = dvdn*ny;
        double divu = dudx + dvdy;
        double txx  = mu_f * (2*dudx - (2.0/3.0)*divu);
        double tyy  = mu_f * (2*dvdy - (2.0/3.0)*divu);
        double txy  = mu_f * (dudy + dvdx);
        double Fv[4] = {
            0.0,
            txx*nx + txy*ny,
            txy*nx + tyy*ny,
            (txx*uf+txy*vf)*nx + (txy*uf+tyy*vf)*ny + kappa_f*dTdn
        };
        if (kL >= 0 && kL < (int)(m.area.size()) && m.area[kL] >= AREA_MIN) {
            double invAL = 1.0/m.area[kL];
            for (int v = 0; v < 4; ++v) R[v*N + kL] += Fv[v]*len2*invAL;
        }
        if (kR >= 0 && kR < (int)(m.area.size()) && m.area[kR] >= AREA_MIN) {
            double invAR = 1.0/m.area[kR];
            for (int v = 0; v < 4; ++v) R[v*N + kR] -= Fv[v]*len2*invAR;
        }
    };

    for (int i = m.i_TEl; i < m.i_TEu; ++i) {
        // Viscous j-fluxes only for airfoil cells (i_TEl..i_TEu-1).
        // Wake cells are inviscid: h₁ there is ~45x smaller (larger outer-boundary
        // separation), making the viscous CFL >> 0.5 on the explicit time step.

        // ---- j=0 wall face: no-slip adiabatic ----
        {
            int kR = i;
            if (m.area[kR] >= AREA_MIN) {
                double nx = m.jfnx[i], ny = m.jfny[i];
                double len2 = m.jflen[i];
                if (len2 >= 1e-14) {
                    double xw = 0.5*(m.node_x(i,0)+m.node_x(i+1,0));
                    double yw = 0.5*(m.node_y(i,0)+m.node_y(i+1,0));
                    double d  = (m.xc[kR]-xw)*nx + (m.yc[kR]-yw)*ny;
                    if (d > 1e-15) {
                        double uc = vel_u(i,0), vc = vel_v(i,0);
                        viscous_face(-1, kR, nx, ny, len2,
                                     uc/d, vc/d, 0.0,   // dTdn=0 (adiabatic)
                                     0.0, 0.0,            // u=v=0 at wall face
                                     mu_eff(i,0), kappa_eff(i,0));
                    }
                }
            }
        }

        // ---- j=1..ncj-1 interior j-faces ----
        for (int j = 1; j < ncj; ++j) {
            int kL = (j-1)*nci+i, kR = j*nci+i;
            if (m.area[kL] < AREA_MIN || m.area[kR] < AREA_MIN) continue;
            double nx = m.jfnx[j*nci+i], ny = m.jfny[j*nci+i];
            double len2 = m.jflen[j*nci+i];
            if (len2 < 1e-14) continue;

            // Distance between cell centers projected onto face normal
            double d_LR = (m.xc[kR]-m.xc[kL])*nx + (m.yc[kR]-m.yc[kL])*ny;
            if (std::abs(d_LR) < 1e-20) continue;

            double uL = vel_u(i,j-1), uR = vel_u(i,j);
            double vL = vel_v(i,j-1), vR = vel_v(i,j);
            double TL = gamma*press(i,j-1)/rho(i,j-1);
            double TR = gamma*press(i,j  )/rho(i,j  );

            double mu_f     = 0.5*(mu_eff(i,j-1)    + mu_eff(i,j));
            double kappa_f  = 0.5*(kappa_eff(i,j-1)  + kappa_eff(i,j));
            viscous_face(kL, kR, nx, ny, len2,
                         (uR-uL)/d_LR, (vR-vL)/d_LR, (TR-TL)/d_LR,
                         0.5*(uL+uR), 0.5*(vL+vR),
                         mu_f, kappa_f);
        }
    }
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

// ---- SA RHS -----------------------------------------------------------------
//
// Computes the residual for the SA equation (ρν̃ transport).
// Thin-layer: j-direction convection + diffusion + source.
// Only for airfoil cells (i_TEl..i_TEu-1); wake cells are inviscid.
//
// Sign: dQ/dt = R_sa  →  R_sa[k] = -(convective) + source + diffusion

// compute_sa_rhs fills:
//   R_sa[k]    = convection + diffusion + production  (no destruction)
//   D_sa[k]    = destruction coefficient (>=0); destruction = -D_sa * SA
//
// The SA update uses a semi-implicit (Jacobi) treatment of the stiff
// destruction term to keep SA >= 0 without tiny time steps:
//   SA_new = max(0, (SA_old + dt * R_sa) / (1 + dt * D_sa))

void Solver::compute_sa_rhs(std::vector<double>& R_sa,
                             std::vector<double>& D_sa,
                             const std::vector<double>& dt_cell) const {
    if (reynolds <= 0.0 || sa_data.empty()) return;

    const Mesh& m = *pmesh;
    std::fill(R_sa.begin(), R_sa.end(), 0.0);
    std::fill(D_sa.begin(), D_sa.end(), 0.0);

    const double mu_mol = 1.0 / reynolds;

    for (int i = m.i_TEl; i < m.i_TEu; ++i) {
        // ---- j-direction: convection + diffusion ----
        for (int j = 0; j <= ncj; ++j) {
            int jL = j-1, jR = j;
            double nx = m.jfnx[j*nci+i], ny = m.jfny[j*nci+i];
            double len = m.jflen[j*nci+i];
            if (len < 1e-14) continue;

            double Q_upwind, v_n_face;

            if (j == 0) {
                // Wall face: SA=0 at wall → ghost = -SA_real.
                // SA value from upstream (upwind): use wall ghost = 0 at face.
                // Convective flux = 0 (v_n ≈ 0 at no-slip wall).
                // Only diffusion through wall (SA=0 gradient drives SA down).
                if (m.area[0*nci+i] < AREA_MIN) continue;
                double xw = 0.5*(m.node_x(i,0)+m.node_x(i+1,0));
                double yw = 0.5*(m.node_y(i,0)+m.node_y(i+1,0));
                double d  = (m.xc[i]-xw)*nx + (m.yc[i]-yw)*ny;
                if (d < 1e-15) continue;
                // Semi-implicit wall drain: SA→0 at wall.  The drain term
                // mu_sa*ρ*ñ/d is proportional to the cell's own SA, so treat it
                // implicitly (add to D_sa) to ensure stability for any dt.
                // Explicit part = mu_sa * ρ * SA_wall / d = 0 (wall SA = 0).
                double nutil = SA(i,0) / rho(i,0);
                double mu_sa = (mu_mol + nutil) / SA_sigma;
                double invA = 1.0 / m.area[i];
                D_sa[i] += mu_sa * rho(i,0) * len * invA / d;
                continue;
            }

            if (j == ncj) continue;  // farfield: SA maintained by BC

            // Interior j-face between cells (i,jL) and (i,jR)
            int kL = jL*nci+i, kR = jR*nci+i;
            if (m.area[kL] < AREA_MIN || m.area[kR] < AREA_MIN) continue;

            // Face normal velocity
            double uL = vel_u(i,jL), vL = vel_v(i,jL);
            double uR = vel_u(i,jR), vR = vel_v(i,jR);
            v_n_face = 0.5*((uL+uR)*nx + (vL+vR)*ny);

            // Upwind SA convective flux: Q_face = Q_upwind
            Q_upwind = (v_n_face >= 0.0) ? SA(i,jL) : SA(i,jR);

            double F_conv = Q_upwind * v_n_face;
            double invAL = 1.0/m.area[kL], invAR = 1.0/m.area[kR];
            R_sa[kL] -= F_conv * len * invAL;
            R_sa[kR] += F_conv * len * invAR;

            // SA diffusion (ν_mol + ν̃)/σ * ∂ν̃/∂n
            double d_LR = (m.xc[kR]-m.xc[kL])*nx + (m.yc[kR]-m.yc[kL])*ny;
            if (std::abs(d_LR) < 1e-8) continue;  // skip TE-region where j-face d_LR≈0
            // Use gradient of Q=ρν̃ (not ν̃) for both diffusion and cb2.
            // With ν̃ = Q/ρ, computing ∂ν̃/∂n from nutil differences generates a spurious
            // term proportional to ∂ρ/∂n even when Q is uniform (e.g. at initialisation).
            // Using ∂Q/∂n instead is zero for uniform SA and recovers the standard
            // (μ+ν̃)/σ * ∂Q/∂n diffusion in the incompressible limit.
            double nutil_f = 0.5*(SA(i,jL)/rho(i,jL) + SA(i,jR)/rho(i,jR));
            double rho_f   = 0.5*(rho(i,jL)+rho(i,jR));
            double d_Qdn   = (SA(i,jR) - SA(i,jL)) / d_LR;
            double mu_sa_f = (mu_mol + nutil_f) / SA_sigma;
            // Semi-implicit diffusion: split into explicit (from neighbour) + implicit (own cell).
            // Avoids explicit diffusion CFL > 0.5 on fine RANS wall cells.
            double diff_coeff = mu_sa_f / d_LR;
            R_sa[kL] += diff_coeff * SA(i,jR) * len * invAL;   // explicit: jR → kL
            R_sa[kR] += diff_coeff * SA(i,jL) * len * invAR;   // explicit: jL → kR
            D_sa[kL] += diff_coeff * len * invAL;               // implicit: kL drains
            D_sa[kR] += diff_coeff * len * invAR;               // implicit: kR drains
            // cb2 volume source: (cb2/σ) * ρ * (∂ν̃/∂n)² ≈ (cb2/σ) * d_Qdn²/ρ per cell.
            // Averaged from adjacent face gradient; split equally to kL and kR.
            // NOTE: no d_LR*invA factor — this is a volume source, not a flux.
            double Fcb2 = SA_cb2 / SA_sigma * d_Qdn * d_Qdn / rho_f;
            R_sa[kL] += 0.5 * Fcb2;
            R_sa[kR] += 0.5 * Fcb2;
        }

        // ---- SA source terms for each airfoil cell ----
        for (int j = 0; j < ncj; ++j) {
            int k = j*nci+i;
            if (m.area[k] < AREA_MIN) continue;

            double d_wall = m.wall_dist[k];
            if (d_wall < 1e-30) continue;   // exactly on wall

            double Q     = SA(i,j);          // ρν̃
            double r_loc = rho(i,j);
            double nutil = Q / r_loc;         // ν̃ = ρν̃/ρ
            double chi   = nutil / (mu_mol / r_loc);  // χ = ν̃/(μ/ρ) = ρν̃*Re/ρ = Q*Re/r
            if (chi < 0.0) { R_sa[k] -= Q / dt_cell[k]; continue; }  // clip negative

            double fv1_v = sa_fv1(chi);
            double fv2_v = 1.0 - chi / (1.0 + chi * fv1_v);

            // Vorticity |ω| (thin-layer: use j-face velocity gradients)
            double S_vort = 0.0;
            {
                // Average of face below (j) and face above (j+1)
                auto dudn_at = [&](int jf) -> double {
                    double nx2 = m.jfnx[jf*nci+i], ny2 = m.jfny[jf*nci+i];
                    double len2 = m.jflen[jf*nci+i];
                    if (len2 < 1e-14) return 0.0;
                    double duL = vel_u(i, jf-1), duR = vel_u(i, jf);
                    double xcL = m.xc[(jf-1)*nci+i], ycL = m.yc[(jf-1)*nci+i];
                    double xcR = m.xc[jf*nci+i],     ycR = m.yc[jf*nci+i];
                    double d2 = (xcR-xcL)*nx2 + (ycR-ycL)*ny2;
                    if (std::abs(d2) < 1e-8) return 0.0;
                    return (duR - duL) / d2;
                };
                double dudn_lo = (j > 0) ? dudn_at(j)   : vel_u(i,0) / std::max(d_wall, 1e-30);
                double dudn_hi = (j+1 < ncj) ? dudn_at(j+1) : dudn_lo;
                S_vort = std::abs(0.5 * (dudn_lo + dudn_hi));
            }

            // Modified vorticity: S̃ = S + ν̃ fv2 / (κ²d²)
            double Stilde = S_vort + nutil * fv2_v / (SA_kappa*SA_kappa * d_wall*d_wall);
            Stilde = std::max(Stilde, 0.3 * S_vort);  // Bradshaw limiter

            // Production: P = cb1 * S̃ * ρν̃
            double P = SA_cb1 * Stilde * Q;

            // Destruction coefficient (semi-implicit: D = -D_coeff * Q)
            // D_coeff = cw1 * fw * (ν̃/d²)  so that D = -D_coeff * Q
            double r_sa    = std::min(nutil / std::max(Stilde*SA_kappa*SA_kappa*d_wall*d_wall, 1e-30), 10.0);
            double D_coeff = SA_cw1 * sa_fw(r_sa) * (nutil / (d_wall*d_wall));

            // R_sa gets production + convection + diffusion (handled above).
            // Destruction is handled semi-implicitly in run():
            //   SA_new = (SA_old + dt*R_sa) / (1 + dt*D_coeff)
            R_sa[k] += P;
            D_sa[k] += std::max(D_coeff, 0.0);
        }
    }
}

// ---- run loop (SSP-RK3 + local time stepping) --------------------------------
// Same formulation as euler2d/src/Solver.cpp but with bridge-cell aware R.

void Solver::run() {
    int N = nci * ncj;
    std::vector<double> R(4*N, 0.0);
    std::vector<double> R_sa(N, 0.0), D_sa(N, 0.0), dt_cell(N, 0.0);
    std::vector<double> Un(4*N, 0.0);
    std::vector<double> SAn(N, 0.0);   // SA backup for RK3

    bool do_sa = (reynolds > 0.0 && !sa_data.empty());

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
    auto pack_sa = [&](std::vector<double>& buf) {
        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                buf[j*nci+i] = SA(i, j);
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
        if (do_sa) pack_sa(SAn);

        // ---- RK stage 1 ----
        apply_bc();
        compute_rhs(R, dt_cell);
        if (do_sa) compute_sa_rhs(R_sa, D_sa, dt_cell);

        double res_rho = 0.0;
        for (int k = 0; k < N; ++k) res_rho += R[k] * R[k];
        res_rho = std::sqrt(res_rho / N);
        if (res0 < 0.0) res0 = res_rho;
        double res_norm = (res0 > 0.0) ? res_rho / res0 : 1.0;

        conv << iter << " " << res_rho << " " << res_norm << "\n";
        if (iter % output_interval == 0 || iter <= 2 || iter == warmup_iters + 1)
            std::cout << "iter " << std::setw(6) << iter
                      << "  res=" << res_rho
                      << "  res/res0=" << res_norm << "\n" << std::flush;

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
                    double u1 = U(vv,i,j), un = Un[vv*N + j*nci+i];
                    U(vv,i,j) = 0.75*un + 0.25*(u1 + dt_cell[j*nci+i]*R[vv*N + j*nci+i]);
                }

        // ---- RK stage 3 ----
        apply_bc();
        compute_rhs(R, dt_cell);

        for (int j = 0; j < ncj; ++j)
            for (int i = 0; i < nci; ++i)
                for (int vv = 0; vv < 4; ++vv) {
                    double u2 = U(vv,i,j), un = Un[vv*N + j*nci+i];
                    U(vv,i,j) = (1.0/3.0)*un + (2.0/3.0)*(u2 + dt_cell[j*nci+i]*R[vv*N + j*nci+i]);
                }

        // ---- SA: forward Euler + semi-implicit destruction (lagged from NS) ----
        // SA is frozen during NS warmup to avoid explosive transient on the initial
        // uniform SA field (SA CFL > 1 at wall cells during stagnation-zone establish).
        // After warmup the NS has established the BL; SA activates from uniform 3/Re.
        if (do_sa && iter > warmup_iters) {
            apply_bc();
            compute_sa_rhs(R_sa, D_sa, dt_cell);
            for (int j = 0; j < ncj; ++j)
                for (int i = 0; i < nci; ++i) {
                    int k = j*nci+i;
                    double dt = dt_cell[k];
                    double denom = 1.0 + dt * D_sa[k];
                    SA(i,j) = std::max(0.0, (SAn[k] + dt*R_sa[k]) / denom);
                }
        }

        if (!std::isfinite(res_rho)) break;

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
