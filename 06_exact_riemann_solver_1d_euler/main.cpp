#include <iostream>

#include "config.hpp"
#include "io.hpp"
#include "solver.hpp"

int main() {
    const PrimitiveState left{
        config::rho_left,
        config::u_left,
        config::p_left
    };

    const PrimitiveState right{
        config::rho_right,
        config::u_right,
        config::p_right
    };

    std::cout << "Exact 1D Euler Riemann solver\n";
    std::cout << "gamma=" << config::gamma_value
              << " nx=" << config::nx
              << " t=" << config::output_time << '\n';
    std::cout << "Left state:  rho=" << left.rho
              << " u=" << left.u
              << " p=" << left.p << '\n';
    std::cout << "Right state: rho=" << right.rho
              << " u=" << right.u
              << " p=" << right.p << '\n';

    const SampledSolution solution =
        compute_exact_solution(left,
                               right,
                               config::gamma_value,
                               config::x_min,
                               config::x_max,
                               config::nx,
                               config::diaphragm_position,
                               config::output_time,
                               config::newton_tolerance,
                               config::newton_max_iterations);

    write_solution_csv("solution.csv", solution);

    std::cout << "Wrote solution.csv\n";
    return 0;
}
