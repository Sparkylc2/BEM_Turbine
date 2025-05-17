
//
// Created by Lukas Campbell on 12/05/2025.
//
#include "../include/headers.h"

namespace BEMSolver {
    using namespace units::literals;


    radian_t newton_raphson(const Blade_Func& f, Blade_Func& df, BladeSection& blade, radian_t phi_0, const dimensionless_t tol, const int max_iter) {

        constexpr radian_t phi_min = 1e-4_rad;
        constexpr radian_t phi_max = radian_t(M_PI_2 - 1e-4);

        radian_t phi = phi_0;
        for (int it = 0; it < max_iter; ++it) {
            const dimensionless_t R  = f(blade, phi);

            if (math::fabs(R) < tol) return phi;

            const dimensionless_t dR = df(blade, phi);

            if (math::fabs(dR) < 1e-12) throw std::runtime_error("dR ≈ 0 in Newton step");

            radian_t dphi { (R / dR).value()};
            dphi = Helpers::clamp(dphi, -0.075_rad, 0.075_rad);

            phi -= dphi;

            if (math::fabs(dphi).value() < tol) return phi;

            phi  = Helpers::clamp(phi, phi_min, phi_max);

        }
        return radian_t(NAN);
    }

    radian_t solve_phi_with_retries(BladeSection& blade, Blade_Func& f, Blade_Func& df) {

        constexpr dimensionless_t k_tol = dimensionless_t(1e-3);
        constexpr int k_max_it = 300;

        constexpr radian_t phi_min = -1*1e-3_rad;
        constexpr radian_t phi_max = radian_t(M_PI_2 - 1e-3);


        const std::vector<radian_t> seeds = {
            math::atan2(dimensionless_t(1.0), blade.g_local_tsr()),
            radian_t(M_PI * 0.25),
            radian_t(M_PI * 0.48),
            phi_min,
            phi_max
        };

        for (const radian_t s : seeds) {
            try {
                if (radian_t phi = newton_raphson(f, df, blade, s, k_tol, k_max_it); !std::isnan(phi.value())) {
                    return phi;
                }
            } catch (...) {/* swallow and try next seed */}
        }

        std::cerr << "solve_phi_with_retries: could not bracket or approximate a root in (0,π/2)" << std::endl;
        std::cerr << "Blade: " << std::endl;
        std::cerr << "Code: " << blade.g_naca_code() << std::endl;
        std::cerr << "Radial Pos: " << blade.g_radial_pos() << std::endl;
        std::cerr << "Twist: " << blade.g_twist_angle() << std::endl;
        std::cerr << "at TSR: " << blade.parent_rotor.g_rotor_tsr() << std::endl;
        return radian_t(NAN);

    }

    scalar_t numerical_derivative(const BladeSection& blade, const radian_t x, const radian_t h) {
        const dimensionless_t f_b = compute_residual(blade, x + h);
        const dimensionless_t f_a = compute_residual(blade, x - h);

        return (f_b - f_a) / (2.0 * h.value());
    }

    void run_bem_solver(BladeSection& blade) {
        // std::cout << "Running BEM solver for blade section at radius: " << blade.g_radial_pos() << "\n";


        Blade_Func f  = [](const BladeSection& b, const radian_t& phi) {
            return compute_residual(b, phi);
        };

        Blade_Func df = [](const BladeSection& b, const radian_t& phi) {
            return numerical_derivative(b, phi, radian_t(1e-8));
        };

        const radian_t phi = solve_phi_with_retries(blade, f, df);
        const auto [a, a_p] = perform_blade_iteration(blade, phi);

        blade.update_induction_factors(a, a_p);

        if (!isnan(a.value()) || !isnan(a_p.value())) {

            std::cout << "BEM solver converged to φ = " << phi << " a: " << a << " a': " << a_p << "\n";
        }

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
        // BEMCorrections::snel_stall_delay_correction(cl, cd, blade);
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
        // BEMCorrections::glauert_buhl_axial_induction(a, F);
        /* ----------------------------------------------------- */

        return {a, a_p};
    }

    dimensionless_t compute_residual(const BladeSection& blade_ref, const radian_t& phi) {

        BladeSection blade = blade_ref;
        const auto [sin_phi, cos_phi] = Helpers::get_clamped_phi_components(phi);
        const auto [a, a_p] = perform_blade_iteration(blade, phi);


        return sin_phi / (1.0 - a) - cos_phi / (blade.g_local_tsr() * (1.0 + a_p));
    }



    void run_bem_backup(BladeSection& blade) {


        /* -------------- initial phi guess ------------- */
        radian_t curr_phi = M_PI_2 * 1_rad - 2.0 / 3.0 * math::atan(1.0 / blade.g_local_tsr());
        /* ---------------------------------------------- */


        /* ---------- calculate phi components ---------- */
        auto [curr_s_phi, curr_c_phi] = Helpers::get_clamped_phi_components(curr_phi);
        /* ---------------------------------------------- */


        /* --- update blade alpha and update cl & cd ---- */
        blade.update_alpha(curr_phi);
        blade.update_cl_cd();
        /* ---------------------------------------------- */


        /* ---------- pull req vars into scope ---------- */
        const dimensionless_t curr_loc_sol = blade.g_local_solidity();
        const dimensionless_t curr_cl = blade.g_cl();
        const dimensionless_t curr_tsr = blade.g_local_tsr();
        /* ---------------------------------------------- */


        /* ---------- initial a and a' guesses ---------- */
        dimensionless_t a = curr_loc_sol * curr_cl * curr_s_phi / (curr_loc_sol * curr_cl * curr_s_phi + 4 * curr_c_phi * curr_c_phi);
        dimensionless_t ap = (1 - 3 * a) / (4 * a - 1);
        /* ---------------------------------------------- */


        /* ------------ enter iteration loop ------------ */
        constexpr size_t max_iter = 5000;
        for (size_t i = 0; i < max_iter; ++i) {


            /* --------------- calc new phi ----------------- */
            curr_phi = math::atan2(dimensionless_t(1.0) - a, curr_tsr * (dimensionless_t(1.0) + ap));
            /* ---------------------------------------------- */


            /* --- update blade alpha and update cl & cd ---- */
            blade.update_alpha(curr_phi);
            blade.update_cl_cd();
            /* ---------------------------------------------- */


            /* -------- calculate new phi components -------- */
            auto [curr_s_phi, curr_c_phi] = Helpers::get_clamped_phi_components(curr_phi);
            /* ---------------------------------------------- */


            /* -------------- get new cl and cd ------------- */
            auto [cl, cd] = std::make_pair(blade.g_cl(), blade.g_cd());
            /* ---------------------------------------------- */


            /* ----------- KARMAN-TSIEN MACH CORRECTION  ----------- */
            BEMCorrections::karman_tsein_mach_correction(blade, cl);
            /* ----------------------------------------------------- */


            /* ------------ SNEL STALL DELAY CORRECTION  ----------- */
            // BEMCorrections::snel_stall_delay_correction(cl, cd, blade);
            /* ----------------------------------------------------- */


            /* --- calc normal and tangential components ---- */
            const dimensionless_t cn = cl * curr_c_phi + cd * curr_s_phi;
            const dimensionless_t ct = cl * curr_s_phi - cd * curr_c_phi;
            /* ---------------------------------------------- */


            /* -------- PRANDTL TIP & ROOT CORRECTION ------- */
            scalar_t F = 1.0;
            BEMCorrections::prandtl_tip_root(F, blade, curr_phi);
            /* ---------------------------------------------- */


            /* --------- calc new a and a' values  ---------- */
            dimensionless_t new_a  = curr_loc_sol * cn / (4.0 * F * curr_s_phi * curr_s_phi + curr_loc_sol * cn);
            dimensionless_t new_ap = curr_loc_sol * ct / (4.0 * F * curr_tsr * curr_s_phi * curr_s_phi) * (1 - new_a);
            /* ---------------------------------------------- */


            /* ------ GLAUERT-BUHL AXIAL INDUCTION CORRECTION ------ */
            // BEMCorrections::glauert_buhl_axial_induction(new_a, F);
            /* ----------------------------------------------------- */


            /* ----------- check if converged ------------- */
            if (std::max(math::abs(new_a - a), math::abs(new_ap - ap)) < 1e-12) break;
            /* ---------------------------------------------- */


            /* ------------ update old a and a' ------------- */
            a = new_a;
            ap = new_ap;
            /* ---------------------------------------------- */
        }

        blade.update_induction_factors(a, ap);

    }


    std::pair<dimensionless_t, dimensionless_t> perform_iteration(const BladeSection& blade_ref, const dimensionless_t old_a, const dimensionless_t old_ap) {


        /* -------------- make local copy --------------- */
        BladeSection blade = blade_ref;
        /* ---------------------------------------------- */


        /* ---------- pull req vars into scope ---------- */
        const dimensionless_t loc_tsr = blade.g_local_tsr();
        const dimensionless_t loc_solidity = blade.g_local_solidity();
        /* ---------------------------------------------- */

        /* --------------- calc new phi ----------------- */
        const radian_t phi = math::atan2(dimensionless_t(1.0) - old_a, loc_tsr * (dimensionless_t(1.0) + old_ap));
        /* ---------------------------------------------- */


        /* --- update blade alpha and update cl & cd ---- */
        blade.update_alpha(phi);
        blade.update_cl_cd();
        /* ---------------------------------------------- */


        /* -------- calculate new phi components -------- */
        auto [curr_s_phi, curr_c_phi] = Helpers::get_clamped_phi_components(phi);
        /* ---------------------------------------------- */


        /* -------------- get new cl and cd ------------- */
        auto [cl, cd] = std::make_pair(blade.g_cl(), blade.g_cd());
        /* ---------------------------------------------- */


        /* --- calc normal and tangential components ---- */
        const dimensionless_t cn = cl * curr_c_phi + cd * curr_s_phi;
        const dimensionless_t ct = cl * curr_s_phi - cd * curr_c_phi;
        /* ---------------------------------------------- */



        /* -------- PRANDTL TIP & ROOT CORRECTION ------- */
        scalar_t F = 1.0;
        BEMCorrections::prandtl_tip_root(F, blade, phi);
        /* ---------------------------------------------- */


        /* --------- calc new a and a' values  ---------- */
        dimensionless_t new_a  = loc_solidity * cn / (loc_solidity * cn + 4 * F * curr_c_phi * curr_c_phi);
        dimensionless_t new_ap = loc_solidity * ct / (4.0 * F * loc_tsr * curr_c_phi * curr_c_phi) * (1 - new_a);
        /* ---------------------------------------------- */

        return {new_a, new_ap};

    }



}
