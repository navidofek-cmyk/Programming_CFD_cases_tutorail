#pragma once

namespace config {

constexpr double H = 1.0;
constexpr double Lx = 32.0 * H;
constexpr double Ly = 8.0 * H;

constexpr double obstacle_x0 = 8.0 * H;
constexpr double obstacle_x1 = 9.0 * H;
constexpr double obstacle_y0 = 3.5 * H;
constexpr double obstacle_y1 = 4.5 * H;

constexpr int Nx = 240;
constexpr int Ny = 120;

constexpr double rho = 1.0;
constexpr double U_in = 1.0;
constexpr double Re = 60.0;
constexpr double nu = U_in * H / Re;

// SIMPLE-like pseudo-time and relaxation parameters.
constexpr int max_iterations = 3000;
constexpr int pressure_iterations = 200;
constexpr double pseudo_dt = 0.02;
constexpr double alpha_u = 0.4;
constexpr double alpha_v = 0.4;
constexpr double alpha_p = 0.15;

constexpr int report_interval = 25;

}  // namespace config
