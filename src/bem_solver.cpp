
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


void BEMSolver::run_bem_solver(BladeSection& blade) {
    std::cout << "Running BEM solver for blade section at radius: " << blade.g_radial_pos() << "\n";


    auto f  = [](const BladeSection& b, const radian_t phi) { return compute_residual(b, phi).value(); };

    auto df = [](const BladeSection& b, const radian_t phi) { return numerical_derivative(b, phi, 1e-6); };

    const radian_t phi = solve_phi_with_retries(blade, f, df);
    const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);

    blade.update_alpha(phi);
    blade.update_cl_cd();

    dimensionless_t c_l = blade.g_cl().value();
    dimensionless_t c_d = blade.g_cd().value();

    /* ----------- KARMAN-TSIEN MACH CORRECTION  ----------- */
    BEMCorrections::karman_tsein_mach_correction(blade, c_l);
    /* ----------------------------------------------------- */


    const dimensionless_t c_n = c_l * cos_phi - c_d * sin_phi;
    const dimensionless_t c_t = c_l * sin_phi + c_d * cos_phi;


    dimensionless_t F = 1.0;

    /* ----------- PRANDTL TIP & ROOT CORRECTION ----------- */
    BEMCorrections::prandtl_tip_root(F, blade, phi);
    /* ----------------------------------------------------- */


    const dimensionless_t sigma_r = blade.g_local_solidity();

    dimensionless_t a = sigma_r * c_n / (4.0 * F * sin_phi * sin_phi + sigma_r * c_n);

    /* ------ GLAUERT-BUHL AXIAL INDUCTION CORRECTION ------ */
    BEMCorrections::glauert_buhl_axial_induction(a, F);
    /* ----------------------------------------------------- */


    const dimensionless_t a_prime = 1.0 / (4.0 * F * sin_phi * cos_phi / (sigma_r * c_t) - 1.0);

    blade.update_induction_factors(dimensionless_t(a),
                                   dimensionless_t(a_prime));

    std::cout << "BEM solver converged to φ = " << phi << " rad\n";
}



std::pair<dimensionless_t, dimensionless_t> BEMSolver::perform_blade_iteration(BladeSection &blade, const radian_t& phi) {

    const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);

    blade.update_alpha(phi);
    blade.update_cl_cd();

    dimensionless_t c_l = blade.g_cl().value();
    dimensionless_t c_d = blade.g_cd().value();
    /* ----------- KARMAN-TSIEN MACH CORRECTION  ----------- */
    BEMCorrections::karman_tsein_mach_correction(blade, c_l);
    /* ----------------------------------------------------- */

    
    dimensionless_t c_n = c_l * cos_phi - c_d * sin_phi;
    dimensionless_t c_t = c_l * sin_phi + c_d * cos_phi;

    scalar_t F = 1.0;
    /* ----------- PRANDTL TIP & ROOT CORRECTION ----------- */
    BEMCorrections::prandtl_tip_root(F, blade, phi);
    /* ----------------------------------------------------- */

    const dimensionless_t local_solidity = blade.g_local_solidity();

    dimensionless_t a   = local_solidity * c_n / (4.0 * F * sin_phi * sin_phi + local_solidity * c_n);
    dimensionless_t a_p = local_solidity * c_t / (4.0 * F * sin_phi * cos_phi - local_solidity * c_t);

    /* ------ GLAUERT-BUHL AXIAL INDUCTION CORRECTION ------ */
    BEMCorrections::glauert_buhl_axial_induction(a, F);
    /* ----------------------------------------------------- */
    return {a, a_p};
}




scalar_t BEMSolver::compute_residual(const BladeSection& blade_ref, const radian_t& phi) {
    
    BladeSection blade = blade_ref;

    const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);
    const auto [a, a_p] = perform_blade_iteration(blade, phi);
    return sin_phi / (1.0 - a) - cos_phi / (blade.g_local_tsr() * (1.0 + a_p));
}
