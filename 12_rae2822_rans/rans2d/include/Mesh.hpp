#pragma once
#include <vector>
#include <string>

struct MeshConfig {
    int    ni_wrap;
    int    nj_normal;
    int    n_airfoil;
    int    n_wake;
    double r_far;
    double tanh_beta;   // j-direction clustering strength (large = finer near wall)
};

struct Mesh {
    int ni, nj;

    std::vector<double> x, y;
    std::vector<double> xc, yc;
    std::vector<double> area;

    std::vector<double> ifnx, ifny, iflen;
    std::vector<double> jfnx, jfny, jflen;

    // Minimum distance from each cell centre to the nearest airfoil wall node.
    // Precomputed in generate(). Size nc_i() * nc_j().
    std::vector<double> wall_dist;

    int i_LE, i_TEl, i_TEu;

    inline double& node_x(int i, int j)       { return x[j * ni + i]; }
    inline double& node_y(int i, int j)       { return y[j * ni + i]; }
    inline const double& node_x(int i, int j) const { return x[j * ni + i]; }
    inline const double& node_y(int i, int j) const { return y[j * ni + i]; }

    inline int nc_i() const { return ni - 1; }
    inline int nc_j() const { return nj - 1; }

    void generate(const MeshConfig& cfg);
    void compute_metrics();
    void write_tecplot(const std::string& filename) const;
};
