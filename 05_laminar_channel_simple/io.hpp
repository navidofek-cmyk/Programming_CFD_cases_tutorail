#pragma once
#include <string>
#include "field.hpp"
#include "mesh.hpp"

void write_csv(const std::string& path, const Mesh& m,
               const Field& u, const Field& v, const Field& p);

// Writes cell-centred u, v, p (and velocity magnitude) to a VTK legacy ASCII
// file readable by ParaView / VisIt.  Points = cell centres.
void write_vtk(const std::string& path, const Mesh& m,
               const Field& u, const Field& v, const Field& p);
