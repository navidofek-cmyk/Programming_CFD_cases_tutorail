#include "IO.hpp"
#include "Solver.hpp"
#include "Mesh.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// ---- mesh VTK export --------------------------------------------------------

void write_mesh_vts(const Mesh& m, const std::string& path) {
    int nci = m.nc_i(), ncj = m.nc_j();
    int ni  = m.ni,     nj  = m.nj;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <StructuredGrid WholeExtent=\"0 " << nci << " 0 " << ncj << " 0 0\">\n"
      << "    <Piece Extent=\"0 " << nci << " 0 " << ncj << " 0 0\">\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni; ++i)
            f << m.node_x(i,j) << " " << m.node_y(i,j) << " 0\n";
    f << "        </DataArray>\n"
      << "      </Points>\n";

    f << "      <CellData Scalars=\"area\">\n";

    // Cell area
    f << "        <DataArray type=\"Float64\" Name=\"area\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            f << m.area[j * nci + i] << "\n";
    f << "        </DataArray>\n";

    // Aspect ratio = max(ds_i, ds_j) / min(ds_i, ds_j)
    f << "        <DataArray type=\"Float64\" Name=\"aspect_ratio\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
            double dxi = m.node_x(i+1,j) - m.node_x(i,j);
            double dyi = m.node_y(i+1,j) - m.node_y(i,j);
            double ds_i = std::sqrt(dxi*dxi + dyi*dyi);
            double dxj = m.node_x(i,j+1) - m.node_x(i,j);
            double dyj = m.node_y(i,j+1) - m.node_y(i,j);
            double ds_j = std::sqrt(dxj*dxj + dyj*dyj);
            double ar = (ds_j > 1e-30 && ds_i > 1e-30)
                        ? std::max(ds_i, ds_j) / std::min(ds_i, ds_j) : 0.0;
            f << ar << "\n";
        }
    }
    f << "        </DataArray>\n";

    f << "      </CellData>\n"
      << "    </Piece>\n"
      << "  </StructuredGrid>\n"
      << "</VTKFile>\n";
}

// ---- VTK StructuredGrid writer (.vts) ---------------------------------------
//
// Layout: CellData, real cells only (i in [0,nci-1], j in [0,ncj-1]).
// Point loop: for j=0..nj-1: for i=0..ni-1: x y 0
// Cell loop:  for j=0..ncj-1: for i=0..nci-1: value
// WholeExtent = "0 nci 0 ncj 0 0"  (node counts: nci+1 x ncj+1)

void write_field(const Solver& s, const std::string& path) {
    const Mesh& m = *s.pmesh;
    int nci = s.nci, ncj = s.ncj;
    int ni  = m.ni,  nj  = m.nj;

    double q_inf2    = s.u_inf*s.u_inf + s.v_inf*s.v_inf;
    double cp_denom  = 0.5 * s.rho_inf * q_inf2;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <StructuredGrid WholeExtent=\"0 " << nci << " 0 " << ncj << " 0 0\">\n"
      << "    <Piece Extent=\"0 " << nci << " 0 " << ncj << " 0 0\">\n";

    // Points (nodes)
    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni; ++i)
            f << m.node_x(i,j) << " " << m.node_y(i,j) << " 0\n";
    f << "        </DataArray>\n"
      << "      </Points>\n";

    // CellData
    f << "      <CellData Scalars=\"density\" Vectors=\"velocity\">\n";

    // density
    f << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            f << s.rho(i,j) << "\n";
    f << "        </DataArray>\n";

    // pressure
    f << "        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            f << s.press(i,j) << "\n";
    f << "        </DataArray>\n";

    // velocity (3-component for VTK)
    f << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            f << s.vel_u(i,j) << " " << s.vel_v(i,j) << " 0\n";
    f << "        </DataArray>\n";

    // mach
    f << "        <DataArray type=\"Float64\" Name=\"mach\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i) {
            double u = s.vel_u(i,j), v = s.vel_v(i,j);
            double c = s.sound(i,j);
            f << std::sqrt(u*u + v*v) / c << "\n";
        }
    f << "        </DataArray>\n";

    // cp
    f << "        <DataArray type=\"Float64\" Name=\"cp\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i) {
            double cp = (cp_denom > 0.0) ? (s.press(i,j) - s.p_inf) / cp_denom : 0.0;
            f << cp << "\n";
        }
    f << "        </DataArray>\n";

    // temperature  T = γ·p/ρ  (non-dim, R=1)
    f << "        <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i)
            f << s.gamma * s.press(i,j) / s.rho(i,j) << "\n";
    f << "        </DataArray>\n";

    // entropy  s = ln(p / ρ^γ)
    f << "        <DataArray type=\"Float64\" Name=\"entropy\" format=\"ascii\">\n";
    for (int j = 0; j < ncj; ++j)
        for (int i = 0; i < nci; ++i) {
            double p = s.press(i,j), r = s.rho(i,j);
            f << std::log(p / std::pow(r, s.gamma)) << "\n";
        }
    f << "        </DataArray>\n";

    f << "      </CellData>\n"
      << "    </Piece>\n"
      << "  </StructuredGrid>\n"
      << "</VTKFile>\n";
}

// ---- VTK PolyData of airfoil surface (.vtp) ---------------------------------
//
// A 1D polyline (cells = line segments) along j=0 from i_TEl to i_TEu.
// PointData: cp, pressure, mach_surface.
// Points are the j=0 WALL NODES (not cell centers) for exact geometry.

void write_surface(const Solver& s, const std::string& path) {
    const Mesh& m = *s.pmesh;
    int i0 = m.i_TEl, i1 = m.i_TEu;  // inclusive node range on airfoil
    int n_pts  = i1 - i0 + 1;         // number of nodes
    int n_cells = n_pts - 1;           // number of line segments

    double q_inf2   = s.u_inf*s.u_inf + s.v_inf*s.v_inf;
    double cp_denom = 0.5 * s.rho_inf * q_inf2;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <PolyData>\n"
      << "    <Piece NumberOfPoints=\"" << n_pts
      << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << n_cells
      << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    // Points (wall nodes at j=0)
    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = i0; i <= i1; ++i)
        f << m.node_x(i, 0) << " " << m.node_y(i, 0) << " 0\n";
    f << "        </DataArray>\n"
      << "      </Points>\n";

    // Lines (connectivity)
    f << "      <Lines>\n"
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int k = 0; k < n_cells; ++k)
        f << k << " " << k + 1 << "\n";
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int k = 1; k <= n_cells; ++k) f << 2*k << "\n";
    f << "        </DataArray>\n"
      << "      </Lines>\n";

    // PointData — average adjacent cell values onto nodes
    // Node i has cells (i-1, 0) and (i, 0); use the average (clamp at ends)
    auto node_prim = [&](int i, double& r, double& u, double& v, double& p) {
        int iL = std::max(i - 1, 0);
        int iR = std::min(i, s.nci - 1);
        r = 0.5 * (s.rho  (iL,0) + s.rho  (iR,0));
        u = 0.5 * (s.vel_u(iL,0) + s.vel_u(iR,0));
        v = 0.5 * (s.vel_v(iL,0) + s.vel_v(iR,0));
        p = 0.5 * (s.press(iL,0) + s.press(iR,0));
    };

    f << "      <PointData Scalars=\"cp\">\n";

    // cp
    f << "        <DataArray type=\"Float64\" Name=\"cp\" format=\"ascii\">\n";
    for (int i = i0; i <= i1; ++i) {
        double r, u, v, p; node_prim(i, r, u, v, p);
        f << ((cp_denom > 0.0) ? (p - s.p_inf) / cp_denom : 0.0) << "\n";
    }
    f << "        </DataArray>\n";

    // pressure
    f << "        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n";
    for (int i = i0; i <= i1; ++i) {
        double r, u, v, p; node_prim(i, r, u, v, p);
        f << p << "\n";
    }
    f << "        </DataArray>\n";

    // mach_surface
    f << "        <DataArray type=\"Float64\" Name=\"mach_surface\" format=\"ascii\">\n";
    for (int i = i0; i <= i1; ++i) {
        double r, u, v, p; node_prim(i, r, u, v, p);
        double c = std::sqrt(s.gamma * p / r);
        f << std::sqrt(u*u + v*v) / c << "\n";
    }
    f << "        </DataArray>\n";

    f << "      </PointData>\n"
      << "    </Piece>\n"
      << "  </PolyData>\n"
      << "</VTKFile>\n";
}

// ---- ASCII Cp along airfoil surface -----------------------------------------

void write_cp(const Solver& s, const std::string& path) {
    const Mesh& m = *s.pmesh;
    double q_inf2   = s.u_inf*s.u_inf + s.v_inf*s.v_inf;
    double cp_denom = 0.5 * s.rho_inf * q_inf2;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);
    f << "# x/c  Cp\n";

    // Airfoil cells: i in [i_TEl, i_TEu-1]  (node i_TEl to i_TEu)
    for (int i = m.i_TEl; i < m.i_TEu; ++i) {
        double xc = m.xc[0 * m.nc_i() + i];  // j=0 cell center x
        double cp = (cp_denom > 0.0) ? (s.press(i, 0) - s.p_inf) / cp_denom : 0.0;
        f << xc << " " << cp << "\n";
    }
}
