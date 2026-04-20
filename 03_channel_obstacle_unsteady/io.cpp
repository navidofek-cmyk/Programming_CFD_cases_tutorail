#include "io.hpp"

#include <fstream>
#include <iomanip>

void write_csv(const std::string& filename, const Mesh& mesh, const Field2D& field) {
    std::ofstream out(filename);
    out << std::setprecision(16);

    for (int j = mesh.ny - 1; j >= 0; --j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << field(i, j);
            if (i + 1 < mesh.nx) {
                out << ',';
            }
        }
        out << '\n';
    }
}

void write_vtk(const std::string& filename,
               const Mesh& mesh,
               const Field2D& u,
               const Field2D& v,
               const Field2D& p) {
    std::ofstream out(filename);
    out << std::setprecision(16);

    out << "# vtk DataFile Version 3.0\n";
    out << "2D channel obstacle unsteady flow\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << mesh.nx + 1 << ' ' << mesh.ny + 1 << " 1\n";

    out << "X_COORDINATES " << mesh.nx + 1 << " double\n";
    for (double x : mesh.xf) {
        out << x << '\n';
    }

    out << "Y_COORDINATES " << mesh.ny + 1 << " double\n";
    for (double y : mesh.yf) {
        out << y << '\n';
    }

    out << "Z_COORDINATES 1 double\n";
    out << "0\n";

    out << "CELL_DATA " << mesh.nx * mesh.ny << '\n';

    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << p(i, j) << '\n';
        }
    }

    out << "VECTORS velocity double\n";
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            out << u(i, j) << ' ' << v(i, j) << " 0\n";
        }
    }
}
