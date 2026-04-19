#include "boundary.hpp"

void apply_velocity_bc(const Mesh& m, Field& u, Field& v, double u_in) {
    // Inlet: u = u_in at inlet face (x=0); interior v-faces at x=0.5dx handled by solver
    for (int j = 0; j < m.Ny; ++j) {
        u(0, j) = u_in;
    }

    // Outlet: zero-gradient for u; v at outlet column handled by v_momentum (aP_eff -= aE)
    for (int j = 0; j < m.Ny; ++j) {
        u(m.Nx, j) = u(m.Nx - 1, j);
    }

    // Bottom and top walls (v = 0 on wall faces)
    for (int i = 0; i < m.Nx; ++i) {
        v(i, 0)    = 0.0;
        v(i, m.Ny) = 0.0;
    }
}

void apply_pressure_bc(const Mesh& m, Field& p) {
    // Outlet: Dirichlet p = 0
    for (int j = 0; j < m.Ny; ++j) {
        p(m.Nx - 1, j) = 0.0;
    }
    // Inlet: zero-gradient
    for (int j = 0; j < m.Ny; ++j) {
        p(0, j) = p(1, j);
    }
    // Top/bottom walls: zero-gradient
    for (int i = 0; i < m.Nx; ++i) {
        p(i, 0)        = p(i, 1);
        p(i, m.Ny - 1) = p(i, m.Ny - 2);
    }
}
