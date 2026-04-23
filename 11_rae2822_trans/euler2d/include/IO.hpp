#pragma once
#include <string>

class Solver;

// Write full flow field as VTK StructuredGrid (.vts) for ParaView
void write_field(const Solver& s, const std::string& path);

// Write airfoil surface as VTK PolyData (.vtp) with Cp, pressure, mach
void write_surface(const Solver& s, const std::string& path);

// Write airfoil Cp to plain ASCII (x/c, Cp)
void write_cp(const Solver& s, const std::string& path);
