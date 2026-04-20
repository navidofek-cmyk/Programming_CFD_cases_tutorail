#include "riemann.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

struct PressureFunctionResult {
    double value = 0.0;
    double derivative = 0.0;
};

PressureFunctionResult pressure_function(double p,
                                         const PrimitiveState& state,
                                         double gamma_value) {
    const double a = sound_speed(state, gamma_value);

    if (p > state.p) {
        // Shock branch.
        const double A = 2.0 / ((gamma_value + 1.0) * state.rho);
        const double B = (gamma_value - 1.0) / (gamma_value + 1.0) * state.p;
        const double sqrt_term = std::sqrt(A / (p + B));
        return {
            (p - state.p) * sqrt_term,
            sqrt_term * (1.0 - 0.5 * (p - state.p) / (p + B))
        };
    }

    // Rarefaction branch.
    const double exponent = (gamma_value - 1.0) / (2.0 * gamma_value);
    const double pressure_ratio = p / state.p;
    const double term = std::pow(pressure_ratio, exponent);
    return {
        2.0 * a / (gamma_value - 1.0) * (term - 1.0),
        (1.0 / (state.rho * a)) * std::pow(pressure_ratio, -(gamma_value + 1.0) / (2.0 * gamma_value))
    };
}

double guess_pressure(const PrimitiveState& left,
                      const PrimitiveState& right,
                      double gamma_value) {
    const double a_l = sound_speed(left, gamma_value);
    const double a_r = sound_speed(right, gamma_value);

    // PVRS estimate used only as a starting guess.
    const double p_pv =
        0.5 * (left.p + right.p)
        - 0.125 * (right.u - left.u) * (left.rho + right.rho) * (a_l + a_r);

    return std::max(1.0e-8, p_pv);
}

PrimitiveState sample_left_rarefaction(const PrimitiveState& left,
                                       const StarRegion& star,
                                       double gamma_value,
                                       double self_similar_speed) {
    const double a_l = sound_speed(left, gamma_value);
    const double head_speed = left.u - a_l;

    if (self_similar_speed <= head_speed) {
        return left;
    }

    const double a_star =
        a_l * std::pow(star.p_star / left.p, (gamma_value - 1.0) / (2.0 * gamma_value));
    const double tail_speed = star.u_star - a_star;

    if (self_similar_speed >= tail_speed) {
        return {
            left.rho * std::pow(star.p_star / left.p, 1.0 / gamma_value),
            star.u_star,
            star.p_star
        };
    }

    const double factor = 2.0 / (gamma_value + 1.0)
                          + (gamma_value - 1.0) / ((gamma_value + 1.0) * a_l)
                                * (left.u - self_similar_speed);

    return {
        left.rho * std::pow(factor, 2.0 / (gamma_value - 1.0)),
        2.0 / (gamma_value + 1.0) *
            (a_l + 0.5 * (gamma_value - 1.0) * left.u + self_similar_speed),
        left.p * std::pow(factor, 2.0 * gamma_value / (gamma_value - 1.0))
    };
}

PrimitiveState sample_right_rarefaction(const PrimitiveState& right,
                                        const StarRegion& star,
                                        double gamma_value,
                                        double self_similar_speed) {
    const double a_r = sound_speed(right, gamma_value);
    const double head_speed = right.u + a_r;

    if (self_similar_speed >= head_speed) {
        return right;
    }

    const double a_star =
        a_r * std::pow(star.p_star / right.p, (gamma_value - 1.0) / (2.0 * gamma_value));
    const double tail_speed = star.u_star + a_star;

    if (self_similar_speed <= tail_speed) {
        return {
            right.rho * std::pow(star.p_star / right.p, 1.0 / gamma_value),
            star.u_star,
            star.p_star
        };
    }

    const double factor = 2.0 / (gamma_value + 1.0)
                          - (gamma_value - 1.0) / ((gamma_value + 1.0) * a_r)
                                * (right.u - self_similar_speed);

    return {
        right.rho * std::pow(factor, 2.0 / (gamma_value - 1.0)),
        2.0 / (gamma_value + 1.0) *
            (-a_r + 0.5 * (gamma_value - 1.0) * right.u + self_similar_speed),
        right.p * std::pow(factor, 2.0 * gamma_value / (gamma_value - 1.0))
    };
}

PrimitiveState sample_left_shock(const PrimitiveState& left,
                                 const StarRegion& star,
                                 double gamma_value,
                                 double self_similar_speed) {
    const double pressure_ratio = star.p_star / left.p;
    const double shock_speed =
        left.u - sound_speed(left, gamma_value)
                     * std::sqrt((gamma_value + 1.0) / (2.0 * gamma_value) * pressure_ratio
                                 + (gamma_value - 1.0) / (2.0 * gamma_value));

    if (self_similar_speed <= shock_speed) {
        return left;
    }

    const double density_ratio =
        (pressure_ratio + (gamma_value - 1.0) / (gamma_value + 1.0))
        / ((gamma_value - 1.0) / (gamma_value + 1.0) * pressure_ratio + 1.0);

    return {
        left.rho * density_ratio,
        star.u_star,
        star.p_star
    };
}

PrimitiveState sample_right_shock(const PrimitiveState& right,
                                  const StarRegion& star,
                                  double gamma_value,
                                  double self_similar_speed) {
    const double pressure_ratio = star.p_star / right.p;
    const double shock_speed =
        right.u + sound_speed(right, gamma_value)
                      * std::sqrt((gamma_value + 1.0) / (2.0 * gamma_value) * pressure_ratio
                                  + (gamma_value - 1.0) / (2.0 * gamma_value));

    if (self_similar_speed >= shock_speed) {
        return right;
    }

    const double density_ratio =
        (pressure_ratio + (gamma_value - 1.0) / (gamma_value + 1.0))
        / ((gamma_value - 1.0) / (gamma_value + 1.0) * pressure_ratio + 1.0);

    return {
        right.rho * density_ratio,
        star.u_star,
        star.p_star
    };
}

}  // namespace

StarRegion solve_star_region(const PrimitiveState& left,
                             const PrimitiveState& right,
                             double gamma_value,
                             double tolerance,
                             int max_iterations) {
    double p = guess_pressure(left, right, gamma_value);

    for (int iter = 0; iter < max_iterations; ++iter) {
        const auto f_l = pressure_function(p, left, gamma_value);
        const auto f_r = pressure_function(p, right, gamma_value);

        const double residual = f_l.value + f_r.value + right.u - left.u;
        const double derivative = f_l.derivative + f_r.derivative;
        const double p_new = std::max(1.0e-8, p - residual / derivative);

        const double relative_change =
            2.0 * std::abs(p_new - p) / (p_new + p);

        p = p_new;

        if (relative_change < tolerance) {
            const auto final_l = pressure_function(p, left, gamma_value);
            const auto final_r = pressure_function(p, right, gamma_value);
            return {
                p,
                0.5 * (left.u + right.u + final_r.value - final_l.value)
            };
        }
    }

    throw std::runtime_error("Exact Riemann solver Newton iteration did not converge.");
}

PrimitiveState sample_exact_solution(const PrimitiveState& left,
                                     const PrimitiveState& right,
                                     const StarRegion& star,
                                     double gamma_value,
                                     double x,
                                     double x0,
                                     double time) {
    if (time <= 0.0) {
        return (x < x0) ? left : right;
    }

    const double s = (x - x0) / time;

    if (s <= star.u_star) {
        if (star.p_star <= left.p) {
            return sample_left_rarefaction(left, star, gamma_value, s);
        }
        return sample_left_shock(left, star, gamma_value, s);
    }

    if (star.p_star <= right.p) {
        return sample_right_rarefaction(right, star, gamma_value, s);
    }
    return sample_right_shock(right, star, gamma_value, s);
}
