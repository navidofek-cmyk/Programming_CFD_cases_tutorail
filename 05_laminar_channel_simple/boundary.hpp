#pragma once
#include "field.hpp"
#include "mesh.hpp"

// Inlet: uniform u = u_in, v = 0
// Outlet: zero-gradient u, zero-gradient v, p = 0 (Dirichlet)
// Top/bottom walls: no-slip (u = 0 via ghost, v = 0 on face)

void apply_velocity_bc(const Mesh& m, Field& u, Field& v, double u_in);
void apply_pressure_bc(const Mesh& m, Field& p);
