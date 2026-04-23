#include "Mesh.hpp"
#include "Solver.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <cmath>

// Simple key=value INI parser (ignores comments starting with #)
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
    if (!f) throw std::runtime_error("Cannot open config: " + path);
    return std::string(std::istreambuf_iterator<char>(f),
                       std::istreambuf_iterator<char>());
}

int main(int argc, char* argv[]) {
    std::string config_path = "config/rae2822.ini";
    if (argc >= 2) config_path = argv[1];

    std::string ini = read_file(config_path);

    // ---- mesh parameters ----
    MeshConfig cfg;
    cfg.ni_wrap    = std::stoi(get_param(ini, "ni_wrap",    "256"));
    cfg.nj_normal  = std::stoi(get_param(ini, "nj_normal",  "65"));
    cfg.n_airfoil  = std::stoi(get_param(ini, "n_airfoil",  "129"));
    cfg.n_wake     = std::stoi(get_param(ini, "n_wake",     "64"));
    cfg.r_far      = std::stod(get_param(ini, "r_far",      "25.0"));

    int expected = 2 * cfg.n_wake + cfg.n_airfoil - 1;
    if (cfg.ni_wrap != expected) {
        std::cerr << "Warning: ni_wrap=" << cfg.ni_wrap
                  << " but 2*n_wake+n_airfoil-1=" << expected
                  << ". Using computed value.\n";
        cfg.ni_wrap = expected;
    }

    // ---- solver parameters ----
    double mach          = std::stod(get_param(ini, "mach",           "0.729"));
    double aoa_deg       = std::stod(get_param(ini, "aoa_deg",        "2.31"));
    double gamma         = std::stod(get_param(ini, "gamma",          "1.4"));
    double cfl           = std::stod(get_param(ini, "cfl",            "1.5"));
    int    max_iter      = std::stoi(get_param(ini, "max_iter",       "20000"));
    double residual_drop = std::stod(get_param(ini, "residual_drop",  "1e-5"));
    int    scheme_order  = std::stoi(get_param(ini, "scheme_order",   "2"));
    int    output_interval = std::stoi(get_param(ini, "output_interval", "500"));

    std::filesystem::create_directories("output");

    // ---- Stage 1: mesh ----
    std::cout << "Generating C-grid: ni=" << cfg.ni_wrap
              << " nj=" << cfg.nj_normal
              << " n_airfoil=" << cfg.n_airfoil
              << " n_wake=" << cfg.n_wake
              << " r_far=" << cfg.r_far << "\n";

    Mesh mesh;
    mesh.generate(cfg);

    std::cout << "i_TEl=" << mesh.i_TEl
              << " i_LE=" << mesh.i_LE
              << " i_TEu=" << mesh.i_TEu << "\n";

    // mesh quality
    int nci = mesh.nc_i(), ncj = mesh.nc_j();
    double amin = 1e99, amax = 0.0;
    int neg_count = 0;
    for (double a : mesh.area) {
        if (a < amin) amin = a;
        if (a > amax) amax = a;
        if (a <= 0.0) ++neg_count;
    }

    double dy1_sum = 0.0; int dy1_n = 0;
    for (int i = mesh.i_TEl; i <= mesh.i_TEu; ++i) {
        double dx = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dy = mesh.node_y(i,1) - mesh.node_y(i,0);
        dy1_sum += std::sqrt(dx*dx + dy*dy);
        ++dy1_n;
    }
    double dy1_avg = dy1_sum / dy1_n;

    double ar_max = 0.0;
    for (int i = 0; i < nci; ++i) {
        double dx_i = mesh.node_x(i+1,0) - mesh.node_x(i,0);
        double dy_i = mesh.node_y(i+1,0) - mesh.node_y(i,0);
        double ds_i = std::sqrt(dx_i*dx_i + dy_i*dy_i);
        double dx_j = mesh.node_x(i,1) - mesh.node_x(i,0);
        double dy_j = mesh.node_y(i,1) - mesh.node_y(i,0);
        double ds_j = std::sqrt(dx_j*dx_j + dy_j*dy_j);
        if (ds_j > 0 && ds_i > 0)
            ar_max = std::max(ar_max, ds_i / ds_j);
    }

    std::cout << "\n=== Mesh statistics ===\n";
    std::cout << "  Cells:              " << nci << " x " << ncj
              << " = " << nci*ncj << "\n";
    std::cout << "  Cell area range:    [" << amin << ", " << amax << "]\n";
    if (neg_count > 0)
        std::cerr << "  WARNING: " << neg_count << " non-positive cell areas!\n";
    else
        std::cout << "  Negative areas:     none\n";
    std::cout << "  First wall cell:    " << dy1_avg << " chords (avg normal height)\n";
    std::cout << "  Max AR at wall:     " << ar_max << "\n";

    mesh.write_tecplot("output/mesh.dat");
    std::cout << "\nWrote output/mesh.dat\n";

    // ---- Stage 2: solver data structures + BCs ----
    std::cout << "\n--- Stage 2: initialising solver ---\n";
    std::cout << "  Mach=" << mach << " AoA=" << aoa_deg << "° gamma=" << gamma << "\n";

    Solver solver;
    solver.init(mesh, mach, aoa_deg, gamma,
                cfl, max_iter, residual_drop,
                scheme_order, output_interval);

    solver.print_bc_diagnostics();

    // ---- Stage 3: run solver ----
    std::cout << "\n--- Stage 3: running solver (order=" << scheme_order << ") ---\n";
    solver.run();

    return 0;
}
