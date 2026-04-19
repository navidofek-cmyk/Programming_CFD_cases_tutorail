#pragma once

namespace config {

// Domain
constexpr double Lx = 6.0;
constexpr double Ly = 1.0;

// Grid (pressure cells)
constexpr int Nx = 120;
constexpr int Ny = 40;

// Fluid
constexpr double rho  = 1.0;
constexpr double U_in = 1.0;
constexpr double Re   = 100.0;
constexpr double nu   = U_in * Ly / Re;

// SIMPLE parameters
constexpr int    max_outer        = 4000;
constexpr int    pressure_sweeps  = 80;
constexpr double sor_omega        = 1.8;
constexpr double pseudo_dt        = 0.05;
constexpr double alpha_u          = 0.7;
constexpr double alpha_v          = 0.7;
constexpr double alpha_p          = 0.3;

// Convergence: stop when RMS change in u between outer iterations < tol_u
constexpr double tol_u = 1.0e-7;

constexpr int report_interval = 50;
constexpr int output_interval = 500;

}  // namespace config
