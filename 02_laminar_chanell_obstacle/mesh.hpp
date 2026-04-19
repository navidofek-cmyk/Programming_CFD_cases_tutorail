#pragma once

#include <vector>

class Mesh {
public:
    Mesh(int nx, int ny, double lx, double ly,
         double obstacle_x0, double obstacle_x1,
         double obstacle_y0, double obstacle_y1);

    [[nodiscard]] int idx(int i, int j) const;
    [[nodiscard]] bool is_inside(int i, int j) const;
    [[nodiscard]] bool is_solid(int i, int j) const;
    [[nodiscard]] bool is_fluid(int i, int j) const;

    int nx = 0;
    int ny = 0;
    double lx = 0.0;
    double ly = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    std::vector<double> xc;
    std::vector<double> yc;
    std::vector<int> solid_mask;
};
