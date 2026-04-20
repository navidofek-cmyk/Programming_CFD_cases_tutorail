#pragma once

#include "field.hpp"
#include "mesh.hpp"

struct SolverSettings {
    double rho = 1.0;
    double nu = 0.0;
    double u_in = 1.0;
    double dt = 0.01;
    int n_steps = 100;
    int poisson_iterations = 100;
    int n_piso_corrections = 2;
    int report_interval = 10;
    int write_interval = 50;
};

struct SolverState {
    Field2D u;
    Field2D v;
    Field2D p;

    explicit SolverState(const Mesh& mesh)
        : u(mesh.nx, mesh.ny, 0.0),
          v(mesh.nx, mesh.ny, 0.0),
          p(mesh.nx, mesh.ny, 0.0) {}
};

struct StepReport {
    double u_change = 0.0;
    double v_change = 0.0;
    double p_corr_residual = 0.0;
    double divergence_norm = 0.0;
    double max_velocity = 0.0;
};

StepReport advance_one_time_step(const Mesh& mesh,
                                 const SolverSettings& settings,
                                 SolverState& state);
