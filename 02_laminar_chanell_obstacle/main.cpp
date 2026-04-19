#include <iostream>

#include "boundary.hpp"
#include "config.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    std::cout << "2D laminar channel with square obstacle\n";
    std::cout << "Nx=" << config::Nx
              << " Ny=" << config::Ny
              << " Re=" << config::Re
              << " nu=" << config::nu << '\n';

    const Mesh mesh(
        config::Nx,
        config::Ny,
        config::Lx,
        config::Ly,
        config::obstacle_x0,
        config::obstacle_x1,
        config::obstacle_y0,
        config::obstacle_y1);

    SolverState state(mesh);
    SolverSettings settings{};
    settings.rho = config::rho;
    settings.nu = config::nu;
    settings.u_in = config::U_in;
    settings.pseudo_dt = config::pseudo_dt;
    settings.alpha_u = config::alpha_u;
    settings.alpha_v = config::alpha_v;
    settings.alpha_p = config::alpha_p;
    settings.max_iterations = config::max_iterations;
    settings.pressure_iterations = config::pressure_iterations;
    settings.report_interval = config::report_interval;

    const SolverReport report = run_simple_solver(mesh, settings, state);

    write_csv("u.csv", mesh, state.u);
    write_csv("v.csv", mesh, state.v);
    write_csv("p.csv", mesh, state.p);

    std::cout << "Finished after " << report.iterations << " iterations\n";
    std::cout << "Final mom_u_res=" << report.u_residual
              << " mom_v_res=" << report.v_residual
              << " p_corr_res=" << report.p_corr_residual
              << " max|U|=" << report.max_velocity << '\n';
    std::cout << "Wrote u.csv, v.csv, p.csv\n";
    return 0;
}
