#include "io.hpp"

#include <fstream>
#include <iomanip>

void write_csv(const std::string& filename, const Mesh& mesh, const Field& field) {
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
