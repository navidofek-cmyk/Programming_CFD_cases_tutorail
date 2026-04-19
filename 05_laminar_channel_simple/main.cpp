#include <iostream>
#include <string>

#include "config.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "solver.hpp"

int main() {
    const Mesh mesh(config::Nx, config::Ny, config::Lx, config::Ly);

    SolverSettings s;
    s.rho             = config::rho;
    s.nu              = config::nu;
    s.u_in            = config::U_in;
    s.pseudo_dt       = config::pseudo_dt;
    s.alpha_u         = config::alpha_u;
    s.alpha_v         = config::alpha_v;
    s.alpha_p         = config::alpha_p;
    s.tol_u           = config::tol_u;
    s.max_outer       = config::max_outer;
    s.pressure_sweeps = config::pressure_sweeps;
    s.sor_omega       = config::sor_omega;
    s.report_interval = config::report_interval;

    SolverState state(mesh);

    std::cout << "Re = " << config::Re
              << "  nu = " << config::nu
              << "  grid " << config::Nx << 'x' << config::Ny << '\n';

    const auto report = run_simple(mesh, s, state);

    write_csv("result.csv", mesh, state.u, state.v, state.p);
    write_vtk("result.vtk", mesh, state.u, state.v, state.p);
    std::cout << "Written result.csv / result.vtk  (iter=" << report.iterations
              << "  div=" << report.div_residual << ")\n";

    // Poiseuille reference: u_max = 1.5 * U_in at x far from inlet
    const double u_max_analytic = 1.5 * config::U_in;
    double u_max_sim = 0.0;
    for (int j = 0; j < mesh.Ny; ++j) {
        const double u_c = 0.5 * (state.u(mesh.Nx - 2, j) + state.u(mesh.Nx - 1, j));
        if (u_c > u_max_sim) u_max_sim = u_c;
    }
    std::cout << "Poiseuille check: u_max(outlet) = " << u_max_sim
              << "  analytic = " << u_max_analytic << '\n';

    return 0;
}
