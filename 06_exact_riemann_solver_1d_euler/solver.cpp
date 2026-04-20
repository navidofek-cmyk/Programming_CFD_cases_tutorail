#include "solver.hpp"

SampledSolution compute_exact_solution(const PrimitiveState& left,
                                       const PrimitiveState& right,
                                       double gamma_value,
                                       double x_min,
                                       double x_max,
                                       int nx,
                                       double diaphragm_position,
                                       double output_time,
                                       double tolerance,
                                       int max_iterations) {
    SampledSolution solution;
    solution.x.resize(static_cast<std::size_t>(nx));
    solution.rho.resize(static_cast<std::size_t>(nx));
    solution.u.resize(static_cast<std::size_t>(nx));
    solution.p.resize(static_cast<std::size_t>(nx));
    solution.e.resize(static_cast<std::size_t>(nx));

    const double dx = (x_max - x_min) / static_cast<double>(nx - 1);
    const StarRegion star = solve_star_region(left, right, gamma_value, tolerance, max_iterations);

    for (int i = 0; i < nx; ++i) {
        const double x = x_min + static_cast<double>(i) * dx;
        const PrimitiveState sample =
            sample_exact_solution(left, right, star, gamma_value, x, diaphragm_position, output_time);

        solution.x[static_cast<std::size_t>(i)] = x;
        solution.rho[static_cast<std::size_t>(i)] = sample.rho;
        solution.u[static_cast<std::size_t>(i)] = sample.u;
        solution.p[static_cast<std::size_t>(i)] = sample.p;
        solution.e[static_cast<std::size_t>(i)] = internal_energy(sample, gamma_value);
    }

    return solution;
}
