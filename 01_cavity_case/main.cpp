#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// USER PARAMETERS
// Change these values first when you want to modify the case.
constexpr int nx = 41;
constexpr int ny = 41;
constexpr double dt = 0.001;
constexpr int nt = 500;
constexpr int nit = 50;
constexpr double rho = 1.0;
constexpr double nu = 0.1;
constexpr double u_lid = 1.0;

int idx(int i, int j, int nx) {
    return j * nx + i;
}

void apply_boundary_conditions(std::vector<double>& u,
                               std::vector<double>& v,
                               std::vector<double>& p,
                               int nx,
                               int ny,
                               double u_lid) {
    for (int i = 0; i < nx; ++i) {
        u[idx(i, 0, nx)] = 0.0;
        v[idx(i, 0, nx)] = 0.0;

        u[idx(i, ny - 1, nx)] = u_lid;
        v[idx(i, ny - 1, nx)] = 0.0;
    }

    for (int j = 0; j < ny; ++j) {
        u[idx(0, j, nx)] = 0.0;
        v[idx(0, j, nx)] = 0.0;

        u[idx(nx - 1, j, nx)] = 0.0;
        v[idx(nx - 1, j, nx)] = 0.0;
    }

    // Pressure boundary conditions:
    // dp/dn = 0 on left, right, and bottom walls.
    // p = 0 on the top wall.
    //
    // This is a simple educational choice for the projection step.
    // The fixed pressure on the top boundary removes the arbitrary constant.
    for (int i = 0; i < nx; ++i) {
        p[idx(i, ny - 1, nx)] = 0.0;
    }

    for (int j = 0; j < ny; ++j) {
        p[idx(0, j, nx)] = p[idx(1, j, nx)];
        p[idx(nx - 1, j, nx)] = p[idx(nx - 2, j, nx)];
    }

    for (int i = 0; i < nx; ++i) {
        p[idx(i, 0, nx)] = p[idx(i, 1, nx)];
    }
}

void apply_velocity_boundary_conditions(std::vector<double>& u,
                                        std::vector<double>& v,
                                        int nx,
                                        int ny,
                                        double u_lid) {
    for (int i = 0; i < nx; ++i) {
        u[idx(i, 0, nx)] = 0.0;
        v[idx(i, 0, nx)] = 0.0;

        u[idx(i, ny - 1, nx)] = u_lid;
        v[idx(i, ny - 1, nx)] = 0.0;
    }

    for (int j = 0; j < ny; ++j) {
        u[idx(0, j, nx)] = 0.0;
        v[idx(0, j, nx)] = 0.0;

        u[idx(nx - 1, j, nx)] = 0.0;
        v[idx(nx - 1, j, nx)] = 0.0;
    }
}

void build_rhs(std::vector<double>& b,
               const std::vector<double>& u,
               const std::vector<double>& v,
               int nx,
               int ny,
               double rho,
               double dt,
               double dx,
               double dy) {
    std::fill(b.begin(), b.end(), 0.0);

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            const double du_dx = (u[idx(i + 1, j, nx)] - u[idx(i - 1, j, nx)]) / (2.0 * dx);
            const double dv_dy = (v[idx(i, j + 1, nx)] - v[idx(i, j - 1, nx)]) / (2.0 * dy);
            const double du_dy = (u[idx(i, j + 1, nx)] - u[idx(i, j - 1, nx)]) / (2.0 * dy);
            const double dv_dx = (v[idx(i + 1, j, nx)] - v[idx(i - 1, j, nx)]) / (2.0 * dx);

            b[idx(i, j, nx)] = rho * ((du_dx + dv_dy) / dt
                                      - du_dx * du_dx
                                      - 2.0 * du_dy * dv_dx
                                      - dv_dy * dv_dy);
        }
    }
}

void solve_pressure_poisson(std::vector<double>& p,
                            const std::vector<double>& b,
                            int nx,
                            int ny,
                            double dx,
                            double dy,
                            int nit) {
    std::vector<double> pn(p.size(), 0.0);
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double denom = 2.0 * (dx2 + dy2);

    for (int it = 0; it < nit; ++it) {
        pn = p;

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                p[idx(i, j, nx)] =
                    ((pn[idx(i + 1, j, nx)] + pn[idx(i - 1, j, nx)]) * dy2
                     + (pn[idx(i, j + 1, nx)] + pn[idx(i, j - 1, nx)]) * dx2
                     - b[idx(i, j, nx)] * dx2 * dy2)
                    / denom;
            }
        }

        for (int i = 0; i < nx; ++i) {
            p[idx(i, ny - 1, nx)] = 0.0;
        }
        for (int j = 0; j < ny; ++j) {
            p[idx(0, j, nx)] = p[idx(1, j, nx)];
            p[idx(nx - 1, j, nx)] = p[idx(nx - 2, j, nx)];
        }
        for (int i = 0; i < nx; ++i) {
            p[idx(i, 0, nx)] = p[idx(i, 1, nx)];
        }
    }
}

void update_velocity(std::vector<double>& u,
                     std::vector<double>& v,
                     const std::vector<double>& p,
                     int nx,
                     int ny,
                     double rho,
                     double nu,
                     double dt,
                     double dx,
                     double dy,
                     double u_lid) {
    std::vector<double> un = u;
    std::vector<double> vn = v;

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            u[idx(i, j, nx)] =
                un[idx(i, j, nx)]
                - dt * un[idx(i, j, nx)] * (un[idx(i, j, nx)] - un[idx(i - 1, j, nx)]) / dx
                - dt * vn[idx(i, j, nx)] * (un[idx(i, j, nx)] - un[idx(i, j - 1, nx)]) / dy
                - dt / (2.0 * rho) * (p[idx(i + 1, j, nx)] - p[idx(i - 1, j, nx)]) / dx
                + nu * dt
                      * ((un[idx(i + 1, j, nx)] - 2.0 * un[idx(i, j, nx)] + un[idx(i - 1, j, nx)]) / dx2
                         + (un[idx(i, j + 1, nx)] - 2.0 * un[idx(i, j, nx)] + un[idx(i, j - 1, nx)]) / dy2);

            v[idx(i, j, nx)] =
                vn[idx(i, j, nx)]
                - dt * un[idx(i, j, nx)] * (vn[idx(i, j, nx)] - vn[idx(i - 1, j, nx)]) / dx
                - dt * vn[idx(i, j, nx)] * (vn[idx(i, j, nx)] - vn[idx(i, j - 1, nx)]) / dy
                - dt / (2.0 * rho) * (p[idx(i, j + 1, nx)] - p[idx(i, j - 1, nx)]) / dy
                + nu * dt
                      * ((vn[idx(i + 1, j, nx)] - 2.0 * vn[idx(i, j, nx)] + vn[idx(i - 1, j, nx)]) / dx2
                         + (vn[idx(i, j + 1, nx)] - 2.0 * vn[idx(i, j, nx)] + vn[idx(i, j - 1, nx)]) / dy2);
        }
    }

    apply_velocity_boundary_conditions(u, v, nx, ny, u_lid);
}

double compute_divergence_norm(const std::vector<double>& u,
                               const std::vector<double>& v,
                               int nx,
                               int ny,
                               double dx,
                               double dy) {
    double sum_sq = 0.0;
    int count = 0;

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            const double div =
                (u[idx(i + 1, j, nx)] - u[idx(i - 1, j, nx)]) / (2.0 * dx)
                + (v[idx(i, j + 1, nx)] - v[idx(i, j - 1, nx)]) / (2.0 * dy);
            sum_sq += div * div;
            ++count;
        }
    }

    return (count > 0) ? std::sqrt(sum_sq / static_cast<double>(count)) : 0.0;
}

void write_csv(const std::string& filename,
               const std::vector<double>& field,
               int nx,
               int ny) {
    std::ofstream out(filename);
    out << std::setprecision(16);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            out << field[idx(i, j, nx)];
            if (i + 1 < nx) {
                out << ',';
            }
        }
        out << '\n';
    }
}

void write_vtk(const std::string& filename,
               const std::vector<double>& u,
               const std::vector<double>& v,
               const std::vector<double>& p,
               int nx,
               int ny,
               double dx,
               double dy) {
    std::ofstream out(filename);
    out << std::setprecision(16);

    out << "# vtk DataFile Version 3.0\n";
    out << "2D lid-driven cavity flow\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << nx << ' ' << ny << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << ' ' << dy << " 1\n";
    out << "POINT_DATA " << nx * ny << '\n';

    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            out << p[idx(i, j, nx)] << '\n';
        }
    }

    out << "VECTORS velocity double\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            out << u[idx(i, j, nx)] << ' '
                << v[idx(i, j, nx)] << ' '
                << 0.0 << '\n';
        }
    }
}

double max_abs(const std::vector<double>& field) {
    double value = 0.0;
    for (double x : field) {
        value = std::max(value, std::abs(x));
    }
    return value;
}

int main() {
    const double lx = 1.0;
    const double ly = 1.0;
    const double dx = lx / static_cast<double>(nx - 1);
    const double dy = ly / static_cast<double>(ny - 1);

    std::vector<double> u(nx * ny, 0.0);
    std::vector<double> v(nx * ny, 0.0);
    std::vector<double> p(nx * ny, 0.0);
    std::vector<double> b(nx * ny, 0.0);

    apply_boundary_conditions(u, v, p, nx, ny, u_lid);

    for (int n = 0; n < nt; ++n) {
        // Step 1: build the right-hand side for the pressure Poisson equation.
        // This uses the current velocity field and approximates how far we are
        // from incompressibility after the explicit momentum update.
        build_rhs(b, u, v, nx, ny, rho, dt, dx, dy);

        // Step 2: solve for pressure that will enforce near-zero divergence.
        solve_pressure_poisson(p, b, nx, ny, dx, dy, nit);

        // Step 3: update velocity with explicit convection, explicit diffusion,
        // and pressure-gradient correction. This is a simplified projection-style
        // educational scheme on a collocated uniform grid.
        update_velocity(u, v, p, nx, ny, rho, nu, dt, dx, dy, u_lid);

        if (n % 25 == 0 || n == nt - 1) {
            const double div_norm = compute_divergence_norm(u, v, nx, ny, dx, dy);
            std::cout << "step " << (n + 1) << "/" << nt
                      << "  max|u|=" << max_abs(u)
                      << "  max|v|=" << max_abs(v)
                      << "  div_l2=" << div_norm << '\n';
        }
    }

    write_csv("u.csv", u, nx, ny);
    write_csv("v.csv", v, nx, ny);
    write_csv("p.csv", p, nx, ny);
    write_vtk("cavity.vtk", u, v, p, nx, ny, dx, dy);

    std::cout << "Wrote u.csv, v.csv, p.csv, cavity.vtk\n";
    return 0;
}
