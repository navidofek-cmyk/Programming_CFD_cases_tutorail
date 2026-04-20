#include "boundary.hpp"

namespace {

double average_fluid_neighbor_pressure(const Mesh& mesh, const Field2D& p, int i, int j) {
    double sum = 0.0;
    int count = 0;

    constexpr int di[4] = {1, -1, 0, 0};
    constexpr int dj[4] = {0, 0, 1, -1};

    for (int n = 0; n < 4; ++n) {
        const int ni = i + di[n];
        const int nj = j + dj[n];
        if (mesh.is_inside(ni, nj) && mesh.is_fluid(ni, nj)) {
            sum += p(ni, nj);
            ++count;
        }
    }

    return (count > 0) ? sum / static_cast<double>(count) : 0.0;
}

}  // namespace

void apply_velocity_boundary_conditions(const Mesh& mesh, Field2D& u, Field2D& v, double u_in) {
    for (int j = 0; j < mesh.ny; ++j) {
        u(0, j) = u_in;
        v(0, j) = 0.0;

        u(mesh.nx - 1, j) = u(mesh.nx - 2, j);
        v(mesh.nx - 1, j) = v(mesh.nx - 2, j);
    }

    for (int i = 0; i < mesh.nx; ++i) {
        u(i, 0) = 0.0;
        v(i, 0) = 0.0;

        u(i, mesh.ny - 1) = 0.0;
        v(i, mesh.ny - 1) = 0.0;
    }

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                u(i, j) = 0.0;
                v(i, j) = 0.0;
            }
        }
    }
}

void apply_pressure_boundary_conditions(const Mesh& mesh, Field2D& p) {
    for (int j = 0; j < mesh.ny; ++j) {
        p(0, j) = p(1, j);
        p(mesh.nx - 1, j) = 0.0;
    }

    for (int i = 0; i < mesh.nx; ++i) {
        p(i, 0) = p(i, 1);
        p(i, mesh.ny - 1) = p(i, mesh.ny - 2);
    }

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_solid(i, j)) {
                p(i, j) = average_fluid_neighbor_pressure(mesh, p, i, j);
            }
        }
    }
}
