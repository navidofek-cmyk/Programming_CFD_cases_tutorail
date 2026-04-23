#include "IO.hpp"
#include "Solver.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>

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
