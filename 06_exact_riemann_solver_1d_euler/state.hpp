#pragma once

#include <cmath>
#include <vector>

struct PrimitiveState {
    double rho = 1.0;
    double u = 0.0;
    double p = 1.0;
};

struct SampledSolution {
    std::vector<double> x;
    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> p;
    std::vector<double> e;
};

inline double sound_speed(const PrimitiveState& state, double gamma_value) {
    return std::sqrt(gamma_value * state.p / state.rho);
}

inline double internal_energy(const PrimitiveState& state, double gamma_value) {
    return state.p / ((gamma_value - 1.0) * state.rho);
}
