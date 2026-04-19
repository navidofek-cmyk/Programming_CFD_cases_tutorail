#pragma once

#include "field.hpp"
#include "mesh.hpp"

struct SolverState {
    Field u;
    Field v;
    Field p;

    explicit SolverState(const Mesh& mesh)
        : u(mesh.nx, mesh.ny, 0.0),
          v(mesh.nx, mesh.ny, 0.0),
          p(mesh.nx, mesh.ny, 0.0) {}
};

struct SolverSettings {
    double rho = 1.0;
    double nu = 0.0;
    double u_in = 1.0;
    double pseudo_dt = 0.01;
    double alpha_u = 0.5;
    double alpha_v = 0.5;
    double alpha_p = 0.3;
    int max_iterations = 1000;
    int pressure_iterations = 100;
    int report_interval = 25;
};

struct SolverReport {
    int iterations = 0;
    double u_residual = 0.0;
    double v_residual = 0.0;
    double p_corr_residual = 0.0;
    double max_velocity = 0.0;
};

SolverReport run_simple_solver(const Mesh& mesh, const SolverSettings& settings, SolverState& state);
