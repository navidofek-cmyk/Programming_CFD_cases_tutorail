#pragma once

#include "riemann.hpp"
#include "state.hpp"

SampledSolution compute_exact_solution(const PrimitiveState& left,
                                       const PrimitiveState& right,
                                       double gamma_value,
                                       double x_min,
                                       double x_max,
                                       int nx,
                                       double diaphragm_position,
                                       double output_time,
                                       double tolerance,
                                       int max_iterations);
