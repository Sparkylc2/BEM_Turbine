//
// Created by Lukas Campbell on 15/05/2025.
//

#pragma once
#include "blade_section.h"
#include "headers.h"
#include "rotor.h"

class BladeSection;

namespace BEMCorrections {

    inline void karman_tsein_mach_correction(BladeSection& blade, dimensionless_t& cl) {
        const auto a = blade.g_a();
        const auto ap = blade.g_a_prime();
        const auto omega = blade.parent_rotor.g_angular_velocity();
        const auto v_inf = blade.parent_rotor.g_wind_speed();
        const auto r = blade.g_radial_pos();
        const auto cs = blade.parent_rotor.speed_of_sound;

        const auto w = math::sqrt(math::pow<2>(v_inf * ( 1 - a)) + math::pow<2>(omega * r * (1 + ap)));
        const dimensionless_t mach = w / cs;


        if(mach <= 0.0 || mach >= 0.85) return; // outside valid area

        const auto beta = math::sqrt(math::max(1.0 - mach * mach, 1e-12));
        cl *= (1 + (mach * mach) / (1 + beta)) / beta;
    }


    inline void du_selig_eggers_stall_delay_correction() {

    }



    inline void glauert_buhl_axial_induction(double& a, const double& F) {
    /*
     * From
     * "A New Empirical Relationship Between Thrust Coefficient and
     * Induction Factor for the Turbulent Windmill State” – NREL/TP-500-36834
    */
        if (a > 0.4) {
            const double CT = 4.0 * F * a * (1.0 - a);
            a = 0.5 * (2.0 + CT * (1.0 - 2.0/3.0) - std::sqrt(std::pow(2.0 + CT * (1.0 - 2.0/3.0), 2.0) + 4.0 * (CT / 3.0 - 1.0)));
        }
    }

    inline void prandtl_tip_root(double& F, const BladeSection& blade, const double& phi) {
        const double B = blade.parent_rotor.g_num_blades();
        const double r = blade.g_radial_pos().value();
        const double r_tip = blade.parent_rotor.g_rotor_radius().value();
        const double r_hub = blade.parent_rotor.g_rotor_hub_radius().value();

        const double abs_phi = std::abs(std::sin(phi));

        const double f_tip = 0.5 * B * (r_tip / r - 1.0) / abs_phi;
        const double f_hub = 0.5 * B * (r_hub / r - 1.0) / abs_phi;

        const double F_tip = 2.0 / M_PI * std::acos(std::exp(-f_tip));
        const double F_hub = 2.0 / M_PI * std::acos(std::exp(-f_hub));
        F = F_tip * F_hub;
    }

}



