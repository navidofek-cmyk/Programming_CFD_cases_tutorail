#pragma once
#include "Mesh.hpp"
#include <vector>
#include <string>

class Solver {
public:
    // ---- grid dimensions (real cells, no ghost) ----
    int nci, ncj;

    // ---- NS conservative variables U[var][j+1][i+1]  ----
    // i in [-1, nci], j in [-1, ncj]  →  stored offset by +1
    // var: 0=ρ, 1=ρu, 2=ρv, 3=ρE
    std::vector<double> data;   // 4*(nci+2)*(ncj+2)

    // ---- SA turbulent variable: ρν̃ (same layout as U, 1 variable) ----
    std::vector<double> sa_data;  // (nci+2)*(ncj+2)

    // ---- freestream (non-dim: ρ∞=1, p∞=1/γ, |u∞|=M∞) ----
    double gamma, mach, aoa_deg;
    double rho_inf, u_inf, v_inf, p_inf, E_inf;

    // ---- viscous / turbulence parameters ----
    double reynolds = 0.0;    // 0 = inviscid Euler
    double prandtl  = 0.72;   // molecular Prandtl
    double prandtl_t = 0.9;   // turbulent Prandtl (for SA energy flux)

    // ---- run config ----
    double cfl;               // warmup CFL (1st order)
    double cfl_impl;          // CFL after order switch (2nd order + MUSCL)
    double residual_drop;
    double omega = 1.0;
    int    max_iter, scheme_order, warmup_iters, output_interval;

    double muscl_factor = 0.0;

    const Mesh* pmesh = nullptr;

    void init(const Mesh& m,
              double mach, double aoa_deg, double gamma,
              double reynolds, double prandtl, double prandtl_t,
              double cfl, double cfl_impl, double omega,
              int max_iter, double residual_drop,
              int scheme_order, int warmup_iters, int output_interval);

    void run();

    // ---- NS accessors (with ghost: i=-1..nci, j=-1..ncj) ----
    inline double& U(int var, int i, int j) {
        return data[var*(ncj+2)*(nci+2) + (j+1)*(nci+2) + (i+1)];
    }
    inline const double& U(int var, int i, int j) const {
        return data[var*(ncj+2)*(nci+2) + (j+1)*(nci+2) + (i+1)];
    }

    // ---- SA accessor (with ghost) ----
    inline double& SA(int i, int j) {
        return sa_data[(j+1)*(nci+2) + (i+1)];
    }
    inline const double& SA(int i, int j) const {
        return sa_data[(j+1)*(nci+2) + (i+1)];
    }

    double rho  (int i, int j) const { return U(0, i, j); }
    double vel_u(int i, int j) const { return U(1, i, j) / U(0, i, j); }
    double vel_v(int i, int j) const { return U(2, i, j) / U(0, i, j); }
    double press(int i, int j) const;
    double sound(int i, int j) const;

    // ν̃ = SA/ρ  (turbulent kinematic viscosity, non-dim)
    double nu_t(int i, int j) const { return SA(i,j) / rho(i,j); }
    // μ_t = ρ * ν̃ * fv1  (turbulent dynamic viscosity)
    double mu_t(int i, int j) const;
    // μ_eff = μ + μ_t  (effective dynamic viscosity for viscous fluxes)
    double mu_eff(int i, int j) const;
    // κ_eff for energy equation: κ_lam + κ_turb
    double kappa_eff(int i, int j) const;

    void apply_bc();
    void print_bc_diagnostics() const;

    // SA constants (Spalart-Allmaras 1992, standard 1-equation model)
    static constexpr double SA_cb1   = 0.1355;
    static constexpr double SA_cb2   = 0.622;
    static constexpr double SA_sigma = 2.0 / 3.0;
    static constexpr double SA_kappa = 0.41;
    static constexpr double SA_cw1   = 3.2390;  // cb1/κ² + (1+cb2)/σ
    static constexpr double SA_cw2   = 0.3;
    static constexpr double SA_cw3   = 2.0;
    static constexpr double SA_cv1   = 7.1;

private:
    void set_cons(int i, int j, double r, double u, double v, double p);
    void compute_rhs(std::vector<double>& R, std::vector<double>& dt_cell) const;
    void compute_sa_rhs(std::vector<double>& R_sa,
                        std::vector<double>& D_sa,
                        const std::vector<double>& dt_cell) const;
    void lusgs(const std::vector<double>& R, std::vector<double>& dU) const;

    double sa_fv1(double chi) const;
    double sa_fw(double r) const;
};
