#pragma once

#include "field.hpp"
#include "mesh.hpp"

void apply_velocity_boundary_conditions(const Mesh& mesh, Field& u, Field& v, double u_in);
void apply_pressure_boundary_conditions(const Mesh& mesh, Field& p);
