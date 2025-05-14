
//
// Created by Lukas Campbell on 12/05/2025.
//
#include "../include/headers.h"



bool BEMSolver::bracket(const std::function<double(BladeSection&, double x)> &f, BladeSection& blade, const double min, const double max, const int n, const bool backwards_search, double& b_min, double& b_max) {
    std::vector<double> x = Helpers::linspace(min, max, n);
    if (backwards_search) std::ranges::reverse(x);

    double f_prev = f(blade, x[0]);
    for (size_t i = 1; i < x.size(); ++i) {
        const double f_next = f(blade, x[i]);

        if (f_prev * f_next < 0) {
            if (backwards_search) {
                b_min = x[i - 1];
                b_max = x[i];
            } else {
                b_min = x[i];
                b_max = x[i - 1];
            }
            return true;
        }
        f_prev = f_next;
    }
    b_min = NAN;
    b_max = NAN;
    return false;
}

double BEMSolver::newton_raphson(
        const std::function<double(BladeSection&, double)>&  f,
        const std::function<double(BladeSection&, double)>&  df,
        BladeSection&                                        blade,
        double                                               x0,
        const double                                         tol,
        const int                                            max_iter)
{
    const double phi_min = 0.01;
    const double phi_max = M_PI_2 - 0.01;

    double phi = std::clamp(x0, phi_min, phi_max);

    for (int it = 0; it < max_iter; ++it)
    {
        const double R  = f(blade, phi);
        if (std::fabs(R) < tol) return phi;

        const double dR = df(blade, phi);
        if (std::fabs(dR) < 1e-12)
            throw std::runtime_error("dR ≈ 0 in Newton step");

        double dphi = R / dR;

        dphi = std::clamp(dphi, -0.20, 0.20);

        phi -= dphi;
        phi  = std::clamp(phi, phi_min, phi_max);


        if (std::fabs(dphi) < tol) return phi;
    }

    return NAN;
}




double solve_phi_with_retries(
        BladeSection& blade,
        const std::function<double(BladeSection&,double)>& f,
        const std::function<double(BladeSection&,double)>& df)
{
    constexpr double kTol   = 1e-6;
    constexpr int    kMaxIt = 60;
    constexpr double lo     = 0.01;
    constexpr double hi     = M_PI_2 - 0.01;

    const double V  = blade.parent_rotor.g_wind_speed().value();
    const double r  = blade.g_radial_pos().value();
    const double Ωr = blade.g_angular_vel().value() * r;
    const double axial_guess = std::atan2(V, Ωr);

    const std::vector<double> seeds = {
        axial_guess,
        0.25 * M_PI,
        0.47 * M_PI,
        lo + 1e-3,
        hi - 1e-3
    };

    for (double s : seeds)
    {
        try {
            double φ = BEMSolver::newton_raphson(f, df, blade, s, kTol, kMaxIt);
            if (!std::isnan(φ)) return φ;
        }
        catch (...) {/* swallow and try next seed */}
    }

    // const double phi_lo = 0.01;
    // const double phi_hi = M_PI - 0.01;
    // const int    Nscan  = 800;
    //
    // double x_prev = phi_lo;
    // double f_prev = f(blade, x_prev);
    //
    // double best_x = x_prev;
    // double best_f = std::fabs(f_prev);
    //
    // for (int i = 1; i <= Nscan; ++i)
    // {
    //     const double x = phi_lo + (phi_hi - phi_lo) * i / Nscan;
    //     const double f_cur = f(blade, x);
    //
    //     if (std::fabs(f_cur) < best_f) { best_f = std::fabs(f_cur); best_x = x; }
    //
    //     x_prev = x;
    //     f_prev = f_cur;
    // }
    //
    // const double kTolLoose = 1e-4;
    // if (best_f < kTolLoose)
    //     return best_x;

    throw std::runtime_error(
        "solve_phi_with_retries: could not bracket or approximate a root in (0,π)");
}




double BEMSolver::numerical_derivative(const BladeSection& blade,
                                       double x,
                                       double               h /* = 1e-6 */)
{
    return (compute_residual(blade, radian_t(x + h)).value() -
            compute_residual(blade, radian_t(x - h)).value()) / (2.0 * h);
}


void BEMSolver::run_bem_solver(BladeSection& blade, bool /*eps_all*/)
{
    std::cout << "Running BEM solver for blade section at radius: "
              << unit_cast<double>(blade.g_radial_pos()) << " m\n";

    const double r     = blade.g_radial_pos().value();

    auto f  = [](BladeSection& b, double x)
              { return compute_residual(b, radian_t(x)).value(); };

    auto df = [](BladeSection& b, double x)
              { return numerical_derivative(b, x, 1e-6); };

    const double phi = solve_phi_with_retries(blade, f, df);
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);

    blade.update_alpha(radian_t(phi));
    blade.update_cl_cd();

    const double cl = blade.g_cl().value();
    const double cd = blade.g_cd().value();

    const double cn = cl * cos_phi - cd * sin_phi;
    const double ct = cl * sin_phi + cd * cos_phi;

    const double B      = blade.parent_rotor.g_num_blades().value();
    const double R      = blade.parent_rotor.g_rotor_radius().value();
    const double R_hub  = blade.parent_rotor.g_rotor_hub_radius().value();
    const double F_tip  = 2.0 / M_PI *
                          std::acos(std::exp(-B * (R - r) /
                                            (2.0 * r * std::fabs(sin_phi))));
    const double F_root = 2.0 / M_PI *
                          std::acos(std::exp(-B * (r - R_hub) /
                                            (2.0 * r * std::fabs(sin_phi))));
    const double F      = std::max(1e-4, F_tip * F_root);

    const double sigma_r = blade.g_local_solidity().value();

    double a = sigma_r * cn /
               (4.0 * F * sin_phi * sin_phi + sigma_r * cn);
    if (a > 0.40)
    {
        const double CT = 4.0 * F * a * (1.0 - a);
        a = 0.5 * (2.0 + CT * (1.0 - 2.0/3.0) -
                   std::sqrt(std::pow(2.0 + CT * (1.0 - 2.0/3.0), 2.0) +
                             4.0 * (CT/3.0 - 1.0)));
    }

    const double a_prime = 1.0 / (4.0 * F * sin_phi * cos_phi / (sigma_r * ct) - 1.0);

    blade.update_induction_factors(dimensionless_t(a),
                                   dimensionless_t(a_prime));

    std::cout << "BEM solver converged to φ = " << phi << " rad\n";
}



dimensionless_t BEMSolver::compute_residual(const BladeSection& blade_ref,
                                            const radian_t&     phi)
{
    BladeSection blade = blade_ref;

    // ---------- basic kinematics ------------------------------------------------
    const double V_inf   = blade.parent_rotor.g_wind_speed().value();
    const double r       = blade.g_radial_pos().value();
    const double Omega   = blade.g_angular_vel().value();
    const double lambda_r = (Omega * r) / V_inf;

    const double sphi = std::sin(phi.value());
    const double cphi = std::cos(phi.value());

    const double s = std::copysign(std::max(std::fabs(sphi), 1e-3), sphi);
    const double c = std::copysign(std::max(std::fabs(cphi), 1e-3), cphi);

    // ---------- aerodynamic coefficients ----------------------------------------
    blade.update_alpha(phi);
    blade.update_cl_cd();

    const double cl = blade.g_cl().value();
    const double cd = blade.g_cd().value();

    const double cn = cl * c - cd * s;
    const double ct = cl * s + cd * c;

    // ---------- prandtl loss factors --------------------------------------------
    const double B     = blade.parent_rotor.g_num_blades().value();
    const double R     = blade.parent_rotor.g_rotor_radius().value();
    const double Rhub  = blade.parent_rotor.g_rotor_hub_radius().value();

    const double Ftip  = 2.0/M_PI * std::acos(std::exp(-B*(R - r)/(2.0*r*std::fabs(s))));
    const double Froot = 2.0/M_PI * std::acos(std::exp(-B*(r - Rhub)/(2.0*r*std::fabs(s))));
    const double F = std::max(1e-4, Ftip * Froot);

    const double sigma_r = blade.g_local_solidity().value();



    // ---------- induction factors (with glauerts) --------------------------------
    double a  = sigma_r * cn / (4.0*F*s*s + sigma_r*cn);
    if (a > 0.40) {
        const double CT = 4.0*F*a*(1.0 - a);
        a = 0.5*(2.0 + CT*(1.0 - 2.0/3.0)
              - std::sqrt(std::abs(std::pow(2.0 + CT*(1.0 - 2.0/3.0),2) + 4.0*(CT/3.0 - 1.0))));
    }
    const double a_p = 1.0 /(4.0 * F * s * c / (sigma_r * ct) - 1.0);


    // ---------- residual ---------------------------------------------------------
    const double res = std::sin(phi.value()) / (1.0 - a) - std::cos(phi.value()) / (lambda_r * (1.0 + a_p));
    return dimensionless_t(res);
}
