#include "mesh.hpp"

Mesh::Mesh(int nx_in, int ny_in, double lx_in, double ly_in,
           double obstacle_x0, double obstacle_x1,
           double obstacle_y0, double obstacle_y1)
    : nx(nx_in),
      ny(ny_in),
      lx(lx_in),
      ly(ly_in),
      dx(lx_in / static_cast<double>(nx_in)),
      dy(ly_in / static_cast<double>(ny_in)),
      xc(static_cast<std::size_t>(nx_in), 0.0),
      yc(static_cast<std::size_t>(ny_in), 0.0),
      xf(static_cast<std::size_t>(nx_in + 1), 0.0),
      yf(static_cast<std::size_t>(ny_in + 1), 0.0),
      solid_mask(static_cast<std::size_t>(nx_in) * static_cast<std::size_t>(ny_in), 0) {
    for (int i = 0; i <= nx; ++i) {
        xf[static_cast<std::size_t>(i)] = static_cast<double>(i) * dx;
    }

    for (int j = 0; j <= ny; ++j) {
        yf[static_cast<std::size_t>(j)] = static_cast<double>(j) * dy;
    }

    for (int i = 0; i < nx; ++i) {
        xc[static_cast<std::size_t>(i)] = (static_cast<double>(i) + 0.5) * dx;
    }

    for (int j = 0; j < ny; ++j) {
        yc[static_cast<std::size_t>(j)] = (static_cast<double>(j) + 0.5) * dy;
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const double x = xc[static_cast<std::size_t>(i)];
            const double y = yc[static_cast<std::size_t>(j)];
            const bool in_obstacle =
                x >= obstacle_x0 && x <= obstacle_x1 &&
                y >= obstacle_y0 && y <= obstacle_y1;
            solid_mask[static_cast<std::size_t>(idx(i, j))] = in_obstacle ? 1 : 0;
        }
    }
}

int Mesh::idx(int i, int j) const {
    return j * nx + i;
}

bool Mesh::is_inside(int i, int j) const {
    return i >= 0 && i < nx && j >= 0 && j < ny;
}

bool Mesh::is_solid(int i, int j) const {
    return solid_mask[static_cast<std::size_t>(idx(i, j))] != 0;
}

bool Mesh::is_fluid(int i, int j) const {
    return !is_solid(i, j);
}
