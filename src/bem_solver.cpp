
//
// Created by Lukas Campbell on 12/05/2025.
//
#include "../include/headers.h"

namespace BEMSolver {
    using namespace units::literals;


    radian_t newton_raphson(const Blade_Func& f, Blade_Func& df, BladeSection& blade, radian_t phi_0, const dimensionless_t tol, const int max_iter) {

        constexpr radian_t phi_min = 0.01_rad;
        constexpr radian_t phi_max = radian_t(M_PI_2 - 0.01);

        radian_t phi = Helpers::clamp(phi_0, phi_min, phi_max);

        std::cout << "PHI_0: " << phi << std::endl;
        for (int it = 0; it < max_iter; ++it) {
            const dimensionless_t R  = f(blade, phi);

            std::cout  << "R: " << R << std::endl;
            if (math::fabs(R) < tol) return phi;

            const dimensionless_t dR = df(blade, phi);

            std::cout << "dR: " << dR << std::endl;
            if (math::fabs(dR) < 1e-12) throw std::runtime_error("dR ≈ 0 in Newton step");

            radian_t dphi { (R / dR).value()};
            dphi = Helpers::clamp(dphi, -0.20_rad, 0.20_rad);

            phi -= dphi;
            phi  = Helpers::clamp(phi, phi_min, phi_max);

            std::cout  << "phi: " << phi << std::endl;
            if (math::fabs(dphi).value() < tol) return phi;
        }
        return radian_t(NAN);
    }

    radian_t solve_phi_with_retries(BladeSection& blade, Blade_Func& f, Blade_Func& df) {

        constexpr dimensionless_t k_tol = dimensionless_t(1e-6);
        constexpr int k_max_it = 60;

        constexpr radian_t lo = 0.02_rad;
        constexpr radian_t hi = radian_t(M_PI_2 - 0.02);

        const std::vector<radian_t> seeds = {
            math::atan2(dimensionless_t(1.0), blade.g_local_tsr()),
            radian_t(M_PI * 0.25),
            radian_t(M_PI * 0.47),
            lo,
            hi
        };

        for (const radian_t s : seeds) {
            try {
                if (radian_t phi = newton_raphson(f, df, blade, s, k_tol, k_max_it); !std::isnan(phi.value())) {
                    return phi;
                }
            } catch (...) {/* swallow and try next seed */}
        }

        throw std::runtime_error("solve_phi_with_retries: could not bracket or approximate a root in (0,π)");
    }




    scalar_t numerical_derivative(const BladeSection& blade, const radian_t x, const radian_t h) {
        const dimensionless_t f_b = compute_residual(blade, x + h);
        const dimensionless_t f_a = compute_residual(blade, x - h);

        return (f_b - f_a) / (2.0 * h.value());
    }


    void run_bem_solver(BladeSection& blade) {
        std::cout << "Running BEM solver for blade section at radius: " << blade.g_radial_pos() << "\n";


        Blade_Func f  = [](const BladeSection& b, const radian_t& phi) {
            return compute_residual(b, phi);
        };

        Blade_Func df = [](const BladeSection& b, const radian_t& phi) {
            return numerical_derivative(b, phi, radian_t(1e-6));
        };

        const radian_t phi = solve_phi_with_retries(blade, f, df);
        const auto [a, a_p] = perform_blade_iteration(blade, phi);

        blade.update_induction_factors(a, a_p);

        std::cout << "BEM solver converged to φ = " << phi << " rad\n";
    }



    std::pair<dimensionless_t, dimensionless_t> perform_blade_iteration(BladeSection &blade, const radian_t& phi) {

        const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);

        blade.update_alpha(phi);
        blade.update_cl_cd();

        dimensionless_t cl = blade.g_cl().value();
        dimensionless_t cd = blade.g_cd().value();
        /* ----------- KARMAN-TSIEN MACH CORRECTION  ----------- */
        BEMCorrections::karman_tsein_mach_correction(blade, cl);
        /* ----------------------------------------------------- */

        /* ------------ SNEL STALL DELAY CORRECTION  ----------- */
        BEMCorrections::snel_stall_delay_correction(cl, cd, blade);
        /* ----------------------------------------------------- */
        


        dimensionless_t cn = cl * cos_phi - cd * sin_phi;
        dimensionless_t ct = cl * sin_phi + cd * cos_phi;

        scalar_t F = 1.0;
        /* ----------- PRANDTL TIP & ROOT CORRECTION ----------- */
        BEMCorrections::prandtl_tip_root(F, blade, phi);
        /* ----------------------------------------------------- */

        const dimensionless_t local_solidity = blade.g_local_solidity();

        dimensionless_t a   = local_solidity * cn / (4.0 * F * sin_phi * sin_phi + local_solidity * cn);
        dimensionless_t a_p = local_solidity * ct / (4.0 * F * sin_phi * cos_phi - local_solidity * ct);

        /* ------ GLAUERT-BUHL AXIAL INDUCTION CORRECTION ------ */
        BEMCorrections::glauert_buhl_axial_induction(a, F);
        /* ----------------------------------------------------- */

        return {a, a_p};
    }




    dimensionless_t compute_residual(const BladeSection& blade_ref, const radian_t& phi) {

        BladeSection blade = blade_ref;
        const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);
        const auto [a, a_p] = perform_blade_iteration(blade, phi);
        
        return sin_phi / (1.0 - a) - cos_phi / (blade.g_local_tsr() * (1.0 + a_p));
    }
}
