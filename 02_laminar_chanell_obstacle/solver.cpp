#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "boundary.hpp"

namespace {

double face_u_e(const Field& u, int i, int j) {
    return 0.5 * (u(i, j) + u(i + 1, j));
}

double face_u_w(const Field& u, int i, int j) {
    return 0.5 * (u(i, j) + u(i - 1, j));
}

double face_v_n(const Field& v, int i, int j) {
    return 0.5 * (v(i, j) + v(i, j + 1));
}

double face_v_s(const Field& v, int i, int j) {
    return 0.5 * (v(i, j) + v(i, j - 1));
}

double fluid_neighbor_value(const Mesh& mesh, const Field& f, int ni, int nj, double fallback) {
    if (!mesh.is_inside(ni, nj) || mesh.is_solid(ni, nj)) {
        return fallback;
    }
    return f(ni, nj);
}

double divergence(const Mesh& mesh, const Field& u, const Field& v, int i, int j) {
    const double ue = fluid_neighbor_value(mesh, u, i + 1, j, u(i, j));
    const double uw = fluid_neighbor_value(mesh, u, i - 1, j, u(i, j));
    const double vn = fluid_neighbor_value(mesh, v, i, j + 1, v(i, j));
    const double vs = fluid_neighbor_value(mesh, v, i, j - 1, v(i, j));

    return (ue - uw) / (2.0 * mesh.dx) + (vn - vs) / (2.0 * mesh.dy);
}

double compute_max_velocity(const Mesh& mesh, const Field& u, const Field& v) {
    double max_vel = 0.0;
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }
            const double mag = std::hypot(u(i, j), v(i, j));
            max_vel = std::max(max_vel, mag);
        }
    }
    return max_vel;
}

double l2_change_over_fluid(const Mesh& mesh, const Field& current, const Field& previous) {
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

double solve_pressure_correction(const Mesh& mesh,
                                 const SolverSettings& settings,
                                 const Field& u_star,
                                 const Field& v_star,
                                 Field& p_corr) {
    p_corr.fill(0.0);

    const double ax = 1.0 / (mesh.dx * mesh.dx);
    const double ay = 1.0 / (mesh.dy * mesh.dy);
    const double ap = 2.0 * (ax + ay);
    double rhs_norm = 0.0;

    for (int iter = 0; iter < settings.pressure_iterations; ++iter) {
        rhs_norm = 0.0;

        for (int j = 1; j < mesh.ny - 1; ++j) {
            for (int i = 1; i < mesh.nx - 1; ++i) {
                if (!mesh.is_fluid(i, j)) {
                    continue;
                }

                const double rhs = (settings.rho / settings.pseudo_dt) * divergence(mesh, u_star, v_star, i, j);

                const double pe = fluid_neighbor_value(mesh, p_corr, i + 1, j, p_corr(i, j));
                const double pw = fluid_neighbor_value(mesh, p_corr, i - 1, j, p_corr(i, j));
                const double pn = fluid_neighbor_value(mesh, p_corr, i, j + 1, p_corr(i, j));
                const double ps = fluid_neighbor_value(mesh, p_corr, i, j - 1, p_corr(i, j));

                const double updated = (ax * (pe + pw) + ay * (pn + ps) - rhs) / ap;
                const double delta = updated - p_corr(i, j);
                p_corr(i, j) += 0.8 * delta;

                rhs_norm += rhs * rhs;
            }
        }

        apply_pressure_boundary_conditions(mesh, p_corr);
    }

    int fluid_cells = 0;
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            if (mesh.is_fluid(i, j)) {
                ++fluid_cells;
            }
        }
    }

    return (fluid_cells > 0) ? std::sqrt(rhs_norm / static_cast<double>(fluid_cells)) : 0.0;
}

void predictor_step(const Mesh& mesh,
                    const SolverSettings& settings,
                    const Field& u_old,
                    const Field& v_old,
                    const Field& p,
                    Field& u_star,
                    Field& v_star) {
    const double volume = mesh.dx * mesh.dy;
    const double de = settings.nu * mesh.dy / mesh.dx;
    const double dw = de;
    const double dn = settings.nu * mesh.dx / mesh.dy;
    const double ds = dn;
    const double ap0 = settings.rho * volume / settings.pseudo_dt;

    u_star = u_old;
    v_star = v_old;

    for (int j = 1; j < mesh.ny - 1; ++j) {
        for (int i = 1; i < mesh.nx - 1; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }

            const double ue_face = face_u_e(u_old, i, j);
            const double uw_face = face_u_w(u_old, i, j);
            const double vn_face = face_v_n(v_old, i, j);
            const double vs_face = face_v_s(v_old, i, j);

            const double Fe = settings.rho * ue_face * mesh.dy;
            const double Fw = settings.rho * uw_face * mesh.dy;
            const double Fn = settings.rho * vn_face * mesh.dx;
            const double Fs = settings.rho * vs_face * mesh.dx;

            const double aE = de + std::max(-Fe, 0.0);
            const double aW = dw + std::max(Fw, 0.0);
            const double aN = dn + std::max(-Fn, 0.0);
            const double aS = ds + std::max(Fs, 0.0);
            const double aP = ap0 + aE + aW + aN + aS;

            const double u_e = fluid_neighbor_value(mesh, u_old, i + 1, j, 0.0);
            const double u_w = fluid_neighbor_value(mesh, u_old, i - 1, j, 0.0);
            const double u_n = fluid_neighbor_value(mesh, u_old, i, j + 1, 0.0);
            const double u_s = fluid_neighbor_value(mesh, u_old, i, j - 1, 0.0);

            const double v_e = fluid_neighbor_value(mesh, v_old, i + 1, j, 0.0);
            const double v_w = fluid_neighbor_value(mesh, v_old, i - 1, j, 0.0);
            const double v_n = fluid_neighbor_value(mesh, v_old, i, j + 1, 0.0);
            const double v_s = fluid_neighbor_value(mesh, v_old, i, j - 1, 0.0);

            const double dp_dx = (fluid_neighbor_value(mesh, p, i + 1, j, p(i, j))
                                  - fluid_neighbor_value(mesh, p, i - 1, j, p(i, j)))
                                 / (2.0 * mesh.dx);
            const double dp_dy = (fluid_neighbor_value(mesh, p, i, j + 1, p(i, j))
                                  - fluid_neighbor_value(mesh, p, i, j - 1, p(i, j)))
                                 / (2.0 * mesh.dy);

            const double u_uncorrected =
                (aE * u_e + aW * u_w + aN * u_n + aS * u_s + ap0 * u_old(i, j) - dp_dx * volume) / aP;

            const double v_uncorrected =
                (aE * v_e + aW * v_w + aN * v_n + aS * v_s + ap0 * v_old(i, j) - dp_dy * volume) / aP;

            u_star(i, j) = settings.alpha_u * u_uncorrected + (1.0 - settings.alpha_u) * u_old(i, j);
            v_star(i, j) = settings.alpha_v * v_uncorrected + (1.0 - settings.alpha_v) * v_old(i, j);
        }
    }

    apply_velocity_boundary_conditions(mesh, u_star, v_star, settings.u_in);
}

void correct_velocity_and_pressure(const Mesh& mesh,
                                   const SolverSettings& settings,
                                   const Field& p_corr,
                                   Field& u,
                                   Field& v,
                                   Field& p) {
    for (int j = 1; j < mesh.ny - 1; ++j) {
        for (int i = 1; i < mesh.nx - 1; ++i) {
            if (!mesh.is_fluid(i, j)) {
                continue;
            }

            const double dpc_dx = (fluid_neighbor_value(mesh, p_corr, i + 1, j, p_corr(i, j))
                                   - fluid_neighbor_value(mesh, p_corr, i - 1, j, p_corr(i, j)))
                                  / (2.0 * mesh.dx);
            const double dpc_dy = (fluid_neighbor_value(mesh, p_corr, i, j + 1, p_corr(i, j))
                                   - fluid_neighbor_value(mesh, p_corr, i, j - 1, p_corr(i, j)))
                                  / (2.0 * mesh.dy);

            u(i, j) -= settings.pseudo_dt / settings.rho * dpc_dx;
            v(i, j) -= settings.pseudo_dt / settings.rho * dpc_dy;
            p(i, j) += settings.alpha_p * p_corr(i, j);
        }
    }

    apply_velocity_boundary_conditions(mesh, u, v, settings.u_in);
    apply_pressure_boundary_conditions(mesh, p);
}

}  // namespace

SolverReport run_simple_solver(const Mesh& mesh, const SolverSettings& settings, SolverState& state) {
    SolverReport report{};
    Field u_old(mesh.nx, mesh.ny, 0.0);
    Field v_old(mesh.nx, mesh.ny, 0.0);
    Field p_corr(mesh.nx, mesh.ny, 0.0);

    apply_velocity_boundary_conditions(mesh, state.u, state.v, settings.u_in);
    apply_pressure_boundary_conditions(mesh, state.p);

    for (int iter = 1; iter <= settings.max_iterations; ++iter) {
        u_old = state.u;
        v_old = state.v;

        predictor_step(mesh, settings, u_old, v_old, state.p, state.u, state.v);
        report.p_corr_residual = solve_pressure_correction(mesh, settings, state.u, state.v, p_corr);
        correct_velocity_and_pressure(mesh, settings, p_corr, state.u, state.v, state.p);

        report.iterations = iter;
        report.u_residual = l2_change_over_fluid(mesh, state.u, u_old);
        report.v_residual = l2_change_over_fluid(mesh, state.v, v_old);
        report.max_velocity = compute_max_velocity(mesh, state.u, state.v);

        if (iter % settings.report_interval == 0 || iter == 1 || iter == settings.max_iterations) {
            std::cout << "iter " << iter
                      << "  mom_u_res=" << report.u_residual
                      << "  mom_v_res=" << report.v_residual
                      << "  p_corr_res=" << report.p_corr_residual
                      << "  max|U|=" << report.max_velocity << '\n';
        }
    }

    return report;
}
