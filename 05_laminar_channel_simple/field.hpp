#pragma once
#include <vector>
#include <algorithm>

class Field {
public:
    int nx, ny;

    Field() : nx(0), ny(0) {}
    Field(int nx, int ny, double init = 0.0)
        : nx(nx), ny(ny), data_(static_cast<std::size_t>(nx * ny), init) {}

    double&       operator()(int i, int j)       { return data_[j * nx + i]; }
    const double& operator()(int i, int j) const { return data_[j * nx + i]; }

    void fill(double v) { std::fill(data_.begin(), data_.end(), v); }

    Field& operator=(const Field&) = default;

private:
    std::vector<double> data_;
};
