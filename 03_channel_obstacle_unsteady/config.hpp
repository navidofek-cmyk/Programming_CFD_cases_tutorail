#pragma once

namespace config {

constexpr double H = 1.0;
constexpr double Lx = 32.0 * H;
constexpr double Ly = 8.0 * H;

constexpr double obstacle_x0 = 8.0 * H;
constexpr double obstacle_x1 = 9.0 * H;
constexpr double obstacle_y0 = 3.5 * H;
constexpr double obstacle_y1 = 4.5 * H;

// A slightly smaller default mesh than the steady case keeps the transient
// example easier to run while still showing the obstacle wake clearly.
constexpr int Nx = 160;
constexpr int Ny = 64;

constexpr double rho = 1.0;
constexpr double U_in = 1.0;
constexpr double Re = 60.0;
constexpr double nu = U_in * H / Re;

constexpr double dt = 0.01;
constexpr int n_steps = 600;
constexpr int poisson_iterations = 120;
constexpr int n_piso_corrections = 2;

constexpr int report_interval = 20;
constexpr int write_interval = 100;

}  // namespace config
