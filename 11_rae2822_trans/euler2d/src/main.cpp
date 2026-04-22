#include "Mesh.hpp"
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
        // strip comment
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = line.substr(0, eq);
        // trim k
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

    MeshConfig cfg;
    cfg.ni_wrap    = std::stoi(get_param(ini, "ni_wrap",    "257"));
    cfg.nj_normal  = std::stoi(get_param(ini, "nj_normal",  "65"));
    cfg.n_airfoil  = std::stoi(get_param(ini, "n_airfoil",  "129"));
    cfg.n_wake     = std::stoi(get_param(ini, "n_wake",     "64"));
    cfg.r_far      = std::stod(get_param(ini, "r_far",      "25.0"));

    // Validate ni_wrap consistency
    int expected = 2 * cfg.n_wake + cfg.n_airfoil - 1;
    if (cfg.ni_wrap != expected) {
        std::cerr << "Warning: ni_wrap=" << cfg.ni_wrap
                  << " but 2*n_wake+n_airfoil-1=" << expected
                  << ". Using computed value.\n";
        cfg.ni_wrap = expected;
    }

    std::filesystem::create_directories("output");

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

    // Sanity check: minimum cell area
    double amin = 1e99, amax = -1e99;
    int neg_count = 0;
    for (double a : mesh.area) {
        if (a < amin) amin = a;
        if (a > amax) amax = a;
        if (a <= 0.0) ++neg_count;
    }
    std::cout << "Cell area range: [" << amin << ", " << amax << "]\n";
    if (neg_count > 0)
        std::cerr << "WARNING: " << neg_count << " non-positive cell areas!\n";
    else
        std::cout << "No negative cell areas.\n";

    mesh.write_tecplot("output/mesh.dat");
    std::cout << "Wrote output/mesh.dat\n";

    return 0;
}
