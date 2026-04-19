#include "io.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

void write_csv(const std::string& path, const Mesh& m,
               const Field& u, const Field& v, const Field& p) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);

    f << std::scientific << std::setprecision(8);
    f << "x,y,u,v,p\n";

    for (int j = 0; j < m.Ny; ++j) {
        const double y = (j + 0.5) * m.dy;
        for (int i = 0; i < m.Nx; ++i) {
            const double x    = (i + 0.5) * m.dx;
            const double u_c  = 0.5 * (u(i, j) + u(i + 1, j));
            const double v_c  = 0.5 * (v(i, j) + v(i, j + 1));
            f << x << ',' << y << ',' << u_c << ',' << v_c << ',' << p(i, j) << '\n';
        }
    }
}

void write_vtk(const std::string& path, const Mesh& m,
               const Field& u, const Field& v, const Field& p) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);

    const int N = m.Nx * m.Ny;

    // VTK legacy ASCII — STRUCTURED_POINTS, points at cell centres
    f << "# vtk DataFile Version 3.0\n"
      << "Channel flow\n"
      << "ASCII\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << m.Nx << " " << m.Ny << " 1\n"
      << "ORIGIN "  << 0.5 * m.dx << " " << 0.5 * m.dy << " 0\n"
      << "SPACING " << m.dx << " " << m.dy << " 1\n"
      << "POINT_DATA " << N << "\n";

    f << std::scientific << std::setprecision(6);

    // Pressure
    f << "SCALARS p double 1\nLOOKUP_TABLE default\n";
    for (int j = 0; j < m.Ny; ++j)
        for (int i = 0; i < m.Nx; ++i)
            f << p(i, j) << "\n";

    // Velocity magnitude
    f << "SCALARS speed double 1\nLOOKUP_TABLE default\n";
    for (int j = 0; j < m.Ny; ++j)
        for (int i = 0; i < m.Nx; ++i) {
            const double uc = 0.5 * (u(i, j) + u(i + 1, j));
            const double vc = 0.5 * (v(i, j) + v(i, j + 1));
            f << std::sqrt(uc * uc + vc * vc) << "\n";
        }

    // Velocity vector (ParaView can show streamlines / glyphs)
    f << "VECTORS velocity double\n";
    for (int j = 0; j < m.Ny; ++j)
        for (int i = 0; i < m.Nx; ++i) {
            const double uc = 0.5 * (u(i, j) + u(i + 1, j));
            const double vc = 0.5 * (v(i, j) + v(i, j + 1));
            f << uc << " " << vc << " 0\n";
        }
}
