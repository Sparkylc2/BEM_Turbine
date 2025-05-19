
//
// Created by Lukas Campbell on 12/05/2025.
//
#include "../include/headers.h"

namespace BEMSolver {
    using namespace units::literals;



    void run_bem_solver(BladeSection& blade) {


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
        constexpr size_t max_iter = 10000;
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
