#include "Mesh.hpp"
#include "Solver.hpp"
#include "IO.hpp"
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
        while (!k.empty() && (k.back()==' ' || k.back()=='\t')) k.pop_back();
        while (!k.empty() && (k.front()==' ' || k.front()=='\t')) k.erase(k.begin());
        if (k != key) continue;
        std::string v = line.substr(eq + 1);
        while (!v.empty() && (v.front()==' ' || v.front()=='\t')) v.erase(v.begin());
        while (!v.empty() && (v.back()==' ' || v.back()=='\t')) v.pop_back();
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

    // ---- mesh params ----
    MeshConfig cfg;
    cfg.ni_wrap    = std::stoi(get_param(ini, "ni_wrap",   "640"));
    cfg.nj_normal  = std::stoi(get_param(ini, "nj_normal", "129"));
    cfg.n_airfoil  = std::stoi(get_param(ini, "n_airfoil", "257"));
    cfg.n_wake     = std::stoi(get_param(ini, "n_wake",    "192"));
    cfg.r_far      = std::stod(get_param(ini, "r_far",     "25.0"));
    cfg.tanh_beta  = std::stod(get_param(ini, "tanh_beta", "7.5"));

    int expected = 2*cfg.n_wake + cfg.n_airfoil - 1;
    if (cfg.ni_wrap != expected) {
        std::cerr << "Warning: ni_wrap=" << cfg.ni_wrap
                  << " expected " << expected << ", using computed.\n";
        cfg.ni_wrap = expected;
    }

    // ---- solver params ----
    double mach          = std::stod(get_param(ini, "mach",           "0.729"));
    double aoa_deg       = std::stod(get_param(ini, "aoa_deg",        "2.31"));
    double gamma         = std::stod(get_param(ini, "gamma",          "1.4"));
    double reynolds      = std::stod(get_param(ini, "reynolds",       "0.0"));
    double prandtl       = std::stod(get_param(ini, "prandtl",        "0.72"));
    double prandtl_t     = std::stod(get_param(ini, "prandtl_t",      "0.9"));
    double cfl           = std::stod(get_param(ini, "cfl",            "1.5"));
    double cfl_impl      = std::stod(get_param(ini, "cfl_impl",       "0.8"));
    double omega         = std::stod(get_param(ini, "omega",          "1.0"));
    int    max_iter      = std::stoi(get_param(ini, "max_iter",       "2000"));
    double residual_drop = std::stod(get_param(ini, "residual_drop",  "1e-6"));
    int    scheme_order  = std::stoi(get_param(ini, "scheme_order",   "2"));
    int    warmup_iters  = std::stoi(get_param(ini, "warmup_iters",   "0"));
    int    output_interval = std::stoi(get_param(ini, "output_interval", "200"));

    std::filesystem::create_directories("output");

    // ---- Stage 1: mesh ----
    std::cout << "Generating RANS C-grid: ni=" << cfg.ni_wrap
              << " nj=" << cfg.nj_normal
              << " n_airfoil=" << cfg.n_airfoil
              << " n_wake=" << cfg.n_wake
              << " tanh_beta=" << cfg.tanh_beta << "\n";

    Mesh mesh;
    mesh.generate(cfg);

    int nci = mesh.nc_i(), ncj = mesh.nc_j();
    std::cout << "i_TEl=" << mesh.i_TEl
              << " i_LE=" << mesh.i_LE
              << " i_TEu=" << mesh.i_TEu
              << " cells=" << nci << "x" << ncj << "\n";

    // Mesh quality: y+ estimate
    double dy1_sum = 0.0; int dy1_n = 0;
    for (int i = mesh.i_TEl; i <= mesh.i_TEu; ++i) {
        double dx = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dy = mesh.node_y(i,1) - mesh.node_y(i,0);
        dy1_sum += std::sqrt(dx*dx + dy*dy);
        ++dy1_n;
    }
    double dy1_avg = dy1_sum / dy1_n;
    double re = std::stod(get_param(ini, "reynolds", "6.5e6"));
    double cf = 0.003;
    double u_tau = std::sqrt(cf / 2.0);
    double yplus = dy1_avg * u_tau * re;
    std::cout << "  h1_avg=" << dy1_avg << " c  y+_est=" << yplus << "\n";

    // ---- Stage 2: LU-SGS Euler solver ----
    std::cout << "\n--- Stage 2: LU-SGS Euler solver ---\n";
    std::cout << "  M=" << mach << " AoA=" << aoa_deg
              << "° Re=" << reynolds << " Pr=" << prandtl
              << " CFL_warmup=" << cfl << " CFL_impl=" << cfl_impl << "\n";

    Solver solver;
    solver.init(mesh, mach, aoa_deg, gamma,
                reynolds, prandtl, prandtl_t,
                cfl, cfl_impl, omega, max_iter, residual_drop,
                scheme_order, warmup_iters, output_interval);

    solver.print_bc_diagnostics();
    std::cout << "\n";
    solver.run();

    return 0;
}
