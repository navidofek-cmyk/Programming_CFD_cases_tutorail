#include "solver.hpp"

#include <algorithm>
#include <cmath>

#include "boundary.hpp"

namespace {

double fluid_value_or_self(const Mesh& mesh, const Field2D& field, int i, int j, int ni, int nj) {
    if (!mesh.is_inside(ni, nj) || mesh.is_solid(ni, nj)) {
        return field(i, j);
    }
    return field(ni, nj);
}

double upwind_value(double face_velocity, double cell_value, double neighbor_value) {
    return (face_velocity >= 0.0) ? cell_value : neighbor_value;
}

double compute_divergence(const Mesh& mesh, const Field2D& u, const Field2D& v, int i, int j) {
    const double ue = fluid_value_or_self(mesh, u, i, j, i + 1, j);
    const double uw = fluid_value_or_self(mesh, u, i, j, i - 1, j);
    const double vn = fluid_value_or_self(mesh, v, i, j, i, j + 1);
    const double vs = fluid_value_or_self(mesh, v, i, j, i, j - 1);

    return (ue - uw) / (2.0 * mesh.dx) + (vn - vs) / (2.0 * mesh.dy);
}

double l2_change_on_fluid(const Mesh& mesh, const Field2D& current, const Field2D& previous) {
    double sum = 0.0;
    int count = 0;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }
            const double diff = current(i, j) - previous(i, j);
            sum += diff * diff;
            ++count;
        }
    }

    return (count > 0) ? std::sqrt(sum / static_cast<double>(count)) : 0.0;
}

double divergence_norm(const Mesh& mesh, const Field2D& u, const Field2D& v) {
    double sum = 0.0;
    int count = 0;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }
            const double div = compute_divergence(mesh, u, v, i, j);
            sum += div * div;
            ++count;
        }
    }

    return (count > 0) ? std::sqrt(sum / static_cast<double>(count)) : 0.0;
}

double max_velocity_magnitude(const Mesh& mesh, const Field2D& u, const Field2D& v) {
    double max_value = 0.0;

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }
            max_value = std::max(max_value, std::hypot(u(i, j), v(i, j)));
        }
    }

    return max_value;
}

double solve_pressure_poisson(const Mesh& mesh,
                              const SolverSettings& settings,
                              const Field2D& u_star,
                              const Field2D& v_star,
                              Field2D& p_increment) {
    Field2D p_old(mesh.nx, mesh.ny, 0.0);
    const double dx2 = mesh.dx * mesh.dx;
    const double dy2 = mesh.dy * mesh.dy;
    const double denom = 2.0 * (dx2 + dy2);
    double rhs_l2 = 0.0;
    int fluid_cells = 0;

    p_increment.fill(0.0);

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_fluid(i, j)) {
                ++fluid_cells;
            }
        }
    }

    for (int iter = 0; iter < settings.poisson_iterations; ++iter) {
        p_old = p_increment;
        rhs_l2 = 0.0;

        for (int j = 1; j < mesh.ny - 1; ++j) {
            for (int i = 1; i < mesh.nx - 1; ++i) {
                if (!mesh.is_fluid(i, j)) {
                    continue;
                }

                const double rhs = (settings.rho / settings.dt) * compute_divergence(mesh, u_star, v_star, i, j);

                const double pe = fluid_value_or_self(mesh, p_old, i, j, i + 1, j);
                const double pw = fluid_value_or_self(mesh, p_old, i, j, i - 1, j);
                const double pn = fluid_value_or_self(mesh, p_old, i, j, i, j + 1);
                const double ps = fluid_value_or_self(mesh, p_old, i, j, i, j - 1);

                p_increment(i, j) =
                    ((pe + pw) * dy2 + (pn + ps) * dx2 - rhs * dx2 * dy2) / denom;

                rhs_l2 += rhs * rhs;
            }
        }

        apply_pressure_boundary_conditions(mesh, p_increment);
    }

    return (fluid_cells > 0) ? std::sqrt(rhs_l2 / static_cast<double>(fluid_cells)) : 0.0;
}

void predictor_step(const Mesh& mesh,
                    const SolverSettings& settings,
                    const SolverState& state,
                    Field2D& u_star,
                    Field2D& v_star) {
    u_star = state.u;
    v_star = state.v;

    const double dx2 = mesh.dx * mesh.dx;
    const double dy2 = mesh.dy * mesh.dy;

    for (int j = 1; j < mesh.ny - 1; ++j) {
        for (int i = 1; i < mesh.nx - 1; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }

            const double u_e = fluid_value_or_self(mesh, state.u, i, j, i + 1, j);
            const double u_w = fluid_value_or_self(mesh, state.u, i, j, i - 1, j);
            const double u_n = fluid_value_or_self(mesh, state.u, i, j, i, j + 1);
            const double u_s = fluid_value_or_self(mesh, state.u, i, j, i, j - 1);

            const double v_e = fluid_value_or_self(mesh, state.v, i, j, i + 1, j);
            const double v_w = fluid_value_or_self(mesh, state.v, i, j, i - 1, j);
            const double v_n = fluid_value_or_self(mesh, state.v, i, j, i, j + 1);
            const double v_s = fluid_value_or_self(mesh, state.v, i, j, i, j - 1);

            const double ue_face = 0.5 * (state.u(i, j) + u_e);
            const double uw_face = 0.5 * (state.u(i, j) + u_w);
            const double vn_face = 0.5 * (state.v(i, j) + v_n);
            const double vs_face = 0.5 * (state.v(i, j) + v_s);

            const double conv_u =
                (ue_face * upwind_value(ue_face, state.u(i, j), u_e)
                 - uw_face * upwind_value(uw_face, u_w, state.u(i, j))) / mesh.dx
                + (vn_face * upwind_value(vn_face, state.u(i, j), u_n)
                   - vs_face * upwind_value(vs_face, u_s, state.u(i, j))) / mesh.dy;

            const double conv_v =
                (ue_face * upwind_value(ue_face, state.v(i, j), v_e)
                 - uw_face * upwind_value(uw_face, v_w, state.v(i, j))) / mesh.dx
                + (vn_face * upwind_value(vn_face, state.v(i, j), v_n)
                   - vs_face * upwind_value(vs_face, v_s, state.v(i, j))) / mesh.dy;

            const double diff_u =
                (u_e - 2.0 * state.u(i, j) + u_w) / dx2
                + (u_n - 2.0 * state.u(i, j) + u_s) / dy2;

            const double diff_v =
                (v_e - 2.0 * state.v(i, j) + v_w) / dx2
                + (v_n - 2.0 * state.v(i, j) + v_s) / dy2;

            const double dp_dx =
                (fluid_value_or_self(mesh, state.p, i, j, i + 1, j)
                 - fluid_value_or_self(mesh, state.p, i, j, i - 1, j))
                / (2.0 * mesh.dx);

            const double dp_dy =
                (fluid_value_or_self(mesh, state.p, i, j, i, j + 1)
                 - fluid_value_or_self(mesh, state.p, i, j, i, j - 1))
                / (2.0 * mesh.dy);

            // Educational approximation:
            // explicit convection + diffusion and the old pressure gradient
            // provide a provisional velocity that is then corrected by PISO.
            u_star(i, j) = state.u(i, j) + settings.dt * (-conv_u - dp_dx / settings.rho + settings.nu * diff_u);
            v_star(i, j) = state.v(i, j) + settings.dt * (-conv_v - dp_dy / settings.rho + settings.nu * diff_v);
        }
    }

    apply_velocity_boundary_conditions(mesh, u_star, v_star, settings.u_in);
}

void apply_pressure_correction(const Mesh& mesh,
                               const SolverSettings& settings,
                               const Field2D& p_increment,
                               Field2D& u,
                               Field2D& v,
                               Field2D& p) {
    for (int j = 1; j < mesh.ny - 1; ++j) {
        for (int i = 1; i < mesh.nx - 1; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }

            const double dpc_dx =
                (fluid_value_or_self(mesh, p_increment, i, j, i + 1, j)
                 - fluid_value_or_self(mesh, p_increment, i, j, i - 1, j))
                / (2.0 * mesh.dx);

            const double dpc_dy =
                (fluid_value_or_self(mesh, p_increment, i, j, i, j + 1)
                 - fluid_value_or_self(mesh, p_increment, i, j, i, j - 1))
                / (2.0 * mesh.dy);

            u(i, j) -= settings.dt / settings.rho * dpc_dx;
            v(i, j) -= settings.dt / settings.rho * dpc_dy;
            p(i, j) += p_increment(i, j);
        }
    }

    apply_velocity_boundary_conditions(mesh, u, v, settings.u_in);
    apply_pressure_boundary_conditions(mesh, p);
}

}  // namespace

StepReport advance_one_time_step(const Mesh& mesh,
                                 const SolverSettings& settings,
                                 SolverState& state) {
    StepReport report{};
    Field2D u_old = state.u;
    Field2D v_old = state.v;
    Field2D u_star(mesh.nx, mesh.ny, 0.0);
    Field2D v_star(mesh.nx, mesh.ny, 0.0);
    Field2D p_increment(mesh.nx, mesh.ny, 0.0);

    predictor_step(mesh, settings, state, u_star, v_star);

    state.u = u_star;
    state.v = v_star;

    for (int corr = 0; corr < settings.n_piso_corrections; ++corr) {
        report.p_corr_residual = solve_pressure_poisson(mesh, settings, state.u, state.v, p_increment);
        apply_pressure_correction(mesh, settings, p_increment, state.u, state.v, state.p);
    }

    report.u_change = l2_change_on_fluid(mesh, state.u, u_old);
    report.v_change = l2_change_on_fluid(mesh, state.v, v_old);
    report.divergence_norm = divergence_norm(mesh, state.u, state.v);
    report.max_velocity = max_velocity_magnitude(mesh, state.u, state.v);
    return report;
}
