#pragma once
#include "field.hpp"
#include "mesh.hpp"

struct SolverState {
    Field u;   // (Nx+1) x Ny
    Field v;   // Nx x (Ny+1)
    Field p;   // Nx x Ny

    explicit SolverState(const Mesh& m)
        : u(m.u_nx(), m.u_ny(), 0.0),
          v(m.v_nx(), m.v_ny(), 0.0),
          p(m.p_nx(), m.p_ny(), 0.0) {}
};

struct SolverSettings {
    double rho             = 1.0;
    double nu              = 0.01;
    double u_in            = 1.0;
    double pseudo_dt       = 0.05;
    double alpha_u         = 0.7;
    double alpha_v         = 0.7;
    double alpha_p         = 0.3;
    double tol_u           = 1.0e-7;  // converged when RMS du < tol_u
    int    max_outer       = 2000;
    int    pressure_sweeps = 80;
    double sor_omega       = 1.4;
    int    report_interval = 50;
};

struct SolverReport {
    int    iterations   = 0;
    double div_residual = 0.0;
    double u_change     = 0.0;
    double v_change     = 0.0;
};

SolverReport run_simple(const Mesh& m, const SolverSettings& s, SolverState& state);
