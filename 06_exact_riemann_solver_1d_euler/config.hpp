#pragma once

namespace config {

constexpr double gamma_value = 1.4;

constexpr double x_min = 0.0;
constexpr double x_max = 1.0;
constexpr int nx = 400;

constexpr double diaphragm_position = 0.5;
constexpr double output_time = 0.20;

// Sod shock tube default data.
constexpr double rho_left = 1.0;
constexpr double u_left = 0.0;
constexpr double p_left = 1.0;

constexpr double rho_right = 0.125;
constexpr double u_right = 0.0;
constexpr double p_right = 0.1;

constexpr double newton_tolerance = 1.0e-8;
constexpr int newton_max_iterations = 100;

}  // namespace config
