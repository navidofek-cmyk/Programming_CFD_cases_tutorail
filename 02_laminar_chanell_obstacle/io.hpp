#pragma once

#include <string>

#include "field.hpp"
#include "mesh.hpp"

void write_csv(const std::string& filename, const Mesh& mesh, const Field& field);
