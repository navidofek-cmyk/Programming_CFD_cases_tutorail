#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "boundary.hpp"
#include "config.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    const Mesh mesh(
        config::Nx,
        config::Ny,
        config::Lx,
        config::Ly,
        config::obstacle_x0,
        config::obstacle_x1,
        config::obstacle_y0,
        config::obstacle_y1);

    SolverSettings settings{};
    settings.rho = config::rho;
    settings.nu = config::nu;
    settings.u_in = config::U_in;
    settings.dt = config::dt;
    settings.n_steps = config::n_steps;
    settings.poisson_iterations = config::poisson_iterations;
    settings.n_piso_corrections = config::n_piso_corrections;
    settings.report_interval = config::report_interval;
    settings.write_interval = config::write_interval;

    SolverState state(mesh);
    apply_velocity_boundary_conditions(mesh, state.u, state.v, settings.u_in);
    apply_pressure_boundary_conditions(mesh, state.p);

    std::filesystem::create_directories("output");

    std::cout << "2D unsteady channel with square obstacle\n";
    std::cout << "Nx=" << mesh.nx
              << " Ny=" << mesh.ny
              << " Re=" << config::Re
              << " dt=" << settings.dt
              << " nSteps=" << settings.n_steps << '\n';

    double time = 0.0;

    for (int step = 1; step <= settings.n_steps; ++step) {
        time += settings.dt;
        const StepReport report = advance_one_time_step(mesh, settings, state);

        if (step % settings.report_interval == 0 || step == 1 || step == settings.n_steps) {
            std::cout << "step " << step
                      << "  time=" << time
                      << "  du_l2=" << report.u_change
                      << "  dv_l2=" << report.v_change
                      << "  p_corr_res=" << report.p_corr_residual
                      << "  div_l2=" << report.divergence_norm
                      << "  max|U|=" << report.max_velocity << '\n';
        }

        if (step % settings.write_interval == 0 || step == settings.n_steps) {
            std::ostringstream name;
            name << "output/flow_" << std::setw(5) << std::setfill('0') << step << ".vtk";
            write_vtk(name.str(), mesh, state.u, state.v, state.p);
        }
    }

    write_csv("u.csv", mesh, state.u);
    write_csv("v.csv", mesh, state.v);
    write_csv("p.csv", mesh, state.p);

    std::cout << "Wrote u.csv, v.csv, p.csv and VTK snapshots in output/\n";
    return 0;
}
