#pragma once
#include "Mesh.hpp"
#include <vector>
#include <string>

class Solver {
public:
    // ---- grid dimensions (real cells, no ghost) ----
    int nci, ncj;

    // ---- conservative variables U[var][j+1][i+1]  ----
    // i in [-1, nci], j in [-1, ncj]  →  stored offset by +1
    // var: 0=ρ, 1=ρu, 2=ρv, 3=ρE
    std::vector<double> data;   // 4*(nci+2)*(ncj+2)

    // ---- freestream (non-dim: ρ∞=1, p∞=1/γ, |u∞|=M∞) ----
    double gamma, mach, aoa_deg;
    double rho_inf, u_inf, v_inf, p_inf, E_inf;

    // ---- run config ----
    double cfl;           // CFL during warmup (1st order phase)
    double cfl_impl;      // CFL after warmup (implicit 2nd order phase)
    double residual_drop;
    double omega = 1.0;         // LU-SGS under-relaxation factor
    int    max_iter, scheme_order, warmup_iters, output_interval;

    // MUSCL blending factor [0=first order, 1=full MUSCL]
    double muscl_factor = 0.0;

    const Mesh* pmesh = nullptr;

    void init(const Mesh& m,
              double mach, double aoa_deg, double gamma,
              double cfl, double cfl_impl, double omega,
              int max_iter, double residual_drop,
              int scheme_order, int warmup_iters, int output_interval);

    void run();

    inline double& U(int var, int i, int j) {
        return data[var*(ncj+2)*(nci+2) + (j+1)*(nci+2) + (i+1)];
    }
    inline const double& U(int var, int i, int j) const {
        return data[var*(ncj+2)*(nci+2) + (j+1)*(nci+2) + (i+1)];
    }

    double rho  (int i, int j) const { return U(0, i, j); }
    double vel_u(int i, int j) const { return U(1, i, j) / U(0, i, j); }
    double vel_v(int i, int j) const { return U(2, i, j) / U(0, i, j); }
    double press(int i, int j) const;
    double sound(int i, int j) const;

    void apply_bc();
    void print_bc_diagnostics() const;

private:
    void set_cons(int i, int j, double r, double u, double v, double p);
    void compute_rhs(std::vector<double>& R, std::vector<double>& dt_cell) const;
    void lusgs(const std::vector<double>& R, std::vector<double>& dU) const;
};
