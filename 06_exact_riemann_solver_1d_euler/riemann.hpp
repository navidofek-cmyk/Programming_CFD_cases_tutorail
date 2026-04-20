#pragma once

#include "state.hpp"

struct StarRegion {
    double p_star = 0.0;
    double u_star = 0.0;
};

StarRegion solve_star_region(const PrimitiveState& left,
                             const PrimitiveState& right,
                             double gamma_value,
                             double tolerance,
                             int max_iterations);

PrimitiveState sample_exact_solution(const PrimitiveState& left,
                                     const PrimitiveState& right,
                                     const StarRegion& star,
                                     double gamma_value,
                                     double x,
                                     double x0,
                                     double time);
