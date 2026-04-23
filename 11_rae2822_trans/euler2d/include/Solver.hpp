#pragma once
#include "Mesh.hpp"
#include <vector>
#include <string>

class Solver {
public:
    // ---- grid dimensions (real cells, no ghost) ----
    int nci, ncj;

    // ---- conservative variables U[var][j+1][i+1]  ----
    // i ranges: -1..nci  →  stored as i+1 = 0..nci+1  (total nci+2)
    // j ranges: -1..ncj  →  stored as j+1 = 0..ncj+1  (total ncj+2)
    // var: 0=ρ, 1=ρu, 2=ρv, 3=ρE
    std::vector<double> data;   // 4 * (nci+2) * (ncj+2)

    // ---- freestream (non-dim: ρ∞=1, p∞=1/γ, |u∞|=M∞) ----
    double gamma, mach, aoa_deg;
    double rho_inf, u_inf, v_inf, p_inf, E_inf;

    // ---- run config ----
    double cfl, residual_drop;
    int    max_iter, scheme_order, warmup_iters, output_interval;

    const Mesh* pmesh = nullptr;

    // Initialize field to freestream + apply first BC pass
    void init(const Mesh& m,
              double mach, double aoa_deg, double gamma,
              double cfl, int max_iter, double residual_drop,
              int scheme_order, int warmup_iters, int output_interval);

    // Run the solver loop (SSP-RK3 + local time stepping)
    void run();

    // Conservative variable accessor (supports ghost: i=-1..nci, j=-1..ncj)
    inline double& U(int var, int i, int j) {
        return data[var * (ncj+2) * (nci+2) + (j+1) * (nci+2) + (i+1)];
    }
    inline const double& U(int var, int i, int j) const {
        return data[var * (ncj+2) * (nci+2) + (j+1) * (nci+2) + (i+1)];
    }

    // Primitive helpers (work on any cell including ghost)
    double rho  (int i, int j) const { return U(0, i, j); }
    double vel_u(int i, int j) const { return U(1, i, j) / U(0, i, j); }
    double vel_v(int i, int j) const { return U(2, i, j) / U(0, i, j); }
    double press(int i, int j) const;
    double sound(int i, int j) const;

    // Apply all BCs (wall, farfield, wake cut, i-exit)
    void apply_bc();

    // Print sample ghost values for Stage 2 checkpoint
    void print_bc_diagnostics() const;

private:
    void set_cons(int i, int j, double r, double u, double v, double p);

    // Compute RHS residual R[var*N + j*nci + i] and local dt for real cells.
    // N = nci * ncj
    void compute_rhs(std::vector<double>& R, std::vector<double>& dt_cell) const;
};
