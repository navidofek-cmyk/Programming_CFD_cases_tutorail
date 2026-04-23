#include "Mesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <cmath>

static std::string get_param(const std::string& ini, const std::string& key,
                              const std::string& def = "") {
    std::istringstream ss(ini);
    std::string line;
    while (std::getline(ss, line)) {
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = line.substr(0, eq);
        while (!k.empty() && (k.back() == ' ' || k.back() == '\t')) k.pop_back();
        while (!k.empty() && (k.front() == ' ' || k.front() == '\t')) k.erase(k.begin());
        if (k != key) continue;
        std::string v = line.substr(eq + 1);
        while (!v.empty() && (v.front() == ' ' || v.front() == '\t')) v.erase(v.begin());
        while (!v.empty() && (v.back() == ' ' || v.back() == '\t')) v.pop_back();
        return v;
    }
    return def;
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    return std::string(std::istreambuf_iterator<char>(f), {});
}

int main(int argc, char* argv[]) {
    std::string config_path = "config/rae2822.ini";
    if (argc >= 2) config_path = argv[1];
    std::string ini = read_file(config_path);

    MeshConfig cfg;
    cfg.ni_wrap    = std::stoi(get_param(ini, "ni_wrap",    "513"));
    cfg.nj_normal  = std::stoi(get_param(ini, "nj_normal",  "129"));
    cfg.n_airfoil  = std::stoi(get_param(ini, "n_airfoil",  "257"));
    cfg.n_wake     = std::stoi(get_param(ini, "n_wake",     "128"));
    cfg.r_far      = std::stod(get_param(ini, "r_far",      "25.0"));
    cfg.tanh_beta  = std::stod(get_param(ini, "tanh_beta",  "7.5"));

    int expected = 2 * cfg.n_wake + cfg.n_airfoil - 1;
    if (cfg.ni_wrap != expected) {
        std::cerr << "Warning: ni_wrap=" << cfg.ni_wrap
                  << " expected " << expected << ", using computed.\n";
        cfg.ni_wrap = expected;
    }

    std::filesystem::create_directories("output");

    std::cout << "Generating RANS C-grid: ni=" << cfg.ni_wrap
              << " nj=" << cfg.nj_normal
              << " n_airfoil=" << cfg.n_airfoil
              << " n_wake=" << cfg.n_wake
              << " r_far=" << cfg.r_far
              << " tanh_beta=" << cfg.tanh_beta << "\n";

    Mesh mesh;
    mesh.generate(cfg);

    std::cout << "i_TEl=" << mesh.i_TEl
              << " i_LE=" << mesh.i_LE
              << " i_TEu=" << mesh.i_TEu << "\n";

    int nci = mesh.nc_i(), ncj = mesh.nc_j();

    // ---- mesh quality statistics ----
    double amin = 1e99, amax = 0.0;
    int neg_count = 0;
    for (double a : mesh.area) {
        if (a < amin) amin = a;
        if (a > amax) amax = a;
        if (a <= 0.0) ++neg_count;
    }

    // First-cell normal height at airfoil wall (averaged)
    double dy1_sum = 0.0; int dy1_n = 0;
    for (int i = mesh.i_TEl; i <= mesh.i_TEu; ++i) {
        double dx = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dy = mesh.node_y(i,1) - mesh.node_y(i,0);
        dy1_sum += std::sqrt(dx*dx + dy*dy);
        ++dy1_n;
    }
    double dy1_avg = dy1_sum / dy1_n;
    double dy1_min = 1e99;
    for (int i = mesh.i_TEl; i <= mesh.i_TEu; ++i) {
        double dx = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dy = mesh.node_y(i,1) - mesh.node_y(i,0);
        dy1_min = std::min(dy1_min, std::sqrt(dx*dx + dy*dy));
    }

    // AR at wall
    double ar_max = 0.0;
    for (int i = 0; i < nci; ++i) {
        double dxi = mesh.node_x(i+1,0) - mesh.node_x(i,0);
        double dyi = mesh.node_y(i+1,0) - mesh.node_y(i,0);
        double ds_i = std::sqrt(dxi*dxi + dyi*dyi);
        double dxj = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dyj = mesh.node_y(i,1) - mesh.node_y(i,0);
        double ds_j = std::sqrt(dxj*dxj + dyj*dyj);
        if (ds_j > 0 && ds_i > 0)
            ar_max = std::max(ar_max, ds_i / ds_j);
    }

    // y+ estimate: y+ = y1 * u_tau / nu = y1 * sqrt(Cf/2) * Re
    // Cf ≈ 0.003 for turbulent flat plate at Re=6.5e6
    double re = 6.5e6, cf = 0.003;
    double u_tau = std::sqrt(cf / 2.0);  // non-dim: u_tau / u_inf
    double nu = 1.0 / re;
    double yplus_avg = dy1_avg * u_tau / nu;
    double yplus_min = dy1_min * u_tau / nu;

    std::cout << "\n=== RANS Mesh statistics ===\n";
    std::cout << "  Cells:              " << nci << " x " << ncj
              << " = " << nci*ncj << "\n";
    std::cout << "  Area range:         [" << amin << ", " << amax << "]\n";
    if (neg_count > 0)
        std::cerr << "  WARNING: " << neg_count << " non-positive areas!\n";
    else
        std::cout << "  Negative areas:     none\n";
    std::cout << "  h1 avg (airfoil):   " << dy1_avg << " c\n";
    std::cout << "  h1 min (airfoil):   " << dy1_min << " c\n";
    std::cout << "  y+ estimate avg:    " << yplus_avg << "\n";
    std::cout << "  y+ estimate min:    " << yplus_min << "\n";
    std::cout << "  AR max at wall:     " << ar_max << "\n";

    mesh.write_tecplot("output/mesh.dat");
    std::cout << "\nWrote output/mesh.dat\n";

    return 0;
}
