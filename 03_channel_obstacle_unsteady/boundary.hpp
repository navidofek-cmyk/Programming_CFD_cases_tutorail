#pragma once

#include "field.hpp"
#include "mesh.hpp"

void apply_velocity_boundary_conditions(const Mesh& mesh, Field2D& u, Field2D& v, double u_in);
void apply_pressure_boundary_conditions(const Mesh& mesh, Field2D& p);
