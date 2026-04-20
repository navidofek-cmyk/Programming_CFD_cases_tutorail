#include "io.hpp"

#include <fstream>
#include <iomanip>

void write_solution_csv(const std::string& filename, const SampledSolution& solution) {
    std::ofstream out(filename);
    out << std::setprecision(16);
    out << "x,rho,u,p,e\n";

    for (std::size_t i = 0; i < solution.x.size(); ++i) {
        out << solution.x[i] << ','
            << solution.rho[i] << ','
            << solution.u[i] << ','
            << solution.p[i] << ','
            << solution.e[i] << '\n';
    }
}
