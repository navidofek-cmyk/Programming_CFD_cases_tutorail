#pragma once
#include <vector>
#include <string>
#include <array>

struct MeshConfig {
    int    ni_wrap;       // nodes around airfoil+wake (i-direction, wraps at wake cut)
    int    nj_normal;     // nodes in normal direction (j=0: wall, j=nj-1: farfield)
    int    n_airfoil;     // nodes on airfoil (both surfaces together, LE to TE + TE to LE)
    int    n_wake;        // nodes in each wake branch
    double r_far;         // farfield radius in chords
};

// Node coordinates of the C-grid.
// Indexing: node(i, j) where i in [0, ni_wrap), j in [0, nj_normal).
// node(0,*) = node(ni_wrap-1, *) — the wake-cut is NOT closed here,
// but the two j=0 branches are geometrically coincident on y=0.
struct Mesh {
    int ni, nj;   // node counts

    // Node positions, row-major [j * ni + i]
    std::vector<double> x, y;

    // Cell-center positions (ni-1 x nj-1 cells)
    std::vector<double> xc, yc;

    // Cell areas [j*(ni-1)+i]
    std::vector<double> area;

    // Face normals and lengths (outward from cell (i,j))
    // i-faces: between (i,j) and (i+1,j) — size ni*(nj-1)
    // j-faces: between (i,j) and (i,j+1) — size (ni-1)*nj
    std::vector<double> ifnx, ifny, iflen;  // i-direction faces
    std::vector<double> jfnx, jfny, jflen;  // j-direction faces

    // Indices on the airfoil surface (j=0 nodes)
    int i_LE;          // leading edge node index in i
    int i_TEl;         // trailing edge lower surface (lower wake start)
    int i_TEu;         // trailing edge upper surface (upper wake start)

    inline double& node_x(int i, int j) { return x[j * ni + i]; }
    inline double& node_y(int i, int j) { return y[j * ni + i]; }
    inline const double& node_x(int i, int j) const { return x[j * ni + i]; }
    inline const double& node_y(int i, int j) const { return y[j * ni + i]; }

    inline int nc_i() const { return ni - 1; }
    inline int nc_j() const { return nj - 1; }

    // Generate mesh from config
    void generate(const MeshConfig& cfg);

    // Compute cell centers, areas, face normals (called after generate)
    void compute_metrics();

    // Write Tecplot ASCII (nodes only)
    void write_tecplot(const std::string& filename) const;
};
