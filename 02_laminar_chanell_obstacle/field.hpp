#pragma once

#include <algorithm>
#include <vector>

class Field {
public:
    Field() = default;

    Field(int nx, int ny, double value = 0.0)
        : nx_(nx), ny_(ny), data_(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), value) {}

    void assign(int nx, int ny, double value = 0.0) {
        nx_ = nx;
        ny_ = ny;
        data_.assign(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), value);
    }

    [[nodiscard]] int nx() const { return nx_; }
    [[nodiscard]] int ny() const { return ny_; }

    [[nodiscard]] double& operator()(int i, int j) {
        return data_[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx_) + static_cast<std::size_t>(i)];
    }

    [[nodiscard]] const double& operator()(int i, int j) const {
        return data_[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx_) + static_cast<std::size_t>(i)];
    }

    void fill(double value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    [[nodiscard]] const std::vector<double>& data() const { return data_; }
    [[nodiscard]] std::vector<double>& data() { return data_; }

private:
    int nx_ = 0;
    int ny_ = 0;
    std::vector<double> data_;
};
