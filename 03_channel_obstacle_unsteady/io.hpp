#pragma once

#include <string>

#include "field.hpp"
#include "mesh.hpp"

void write_csv(const std::string& filename, const Mesh& mesh, const Field2D& field);
void write_vtk(const std::string& filename,
               const Mesh& mesh,
               const Field2D& u,
               const Field2D& v,
               const Field2D& p);
