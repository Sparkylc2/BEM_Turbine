//
// Created by Lukas Campbell on 15/05/2025.
//

#pragma once
#include "blade_section.h"
#include "headers.h"
#include "rotor.h"

class BladeSection;

namespace BEMCorrections {
    using namespace units::dimensionless;

    inline void karman_tsein_mach_correction(const BladeSection& blade, dimensionless_t& cl) {
        const auto a = blade.g_a();
        const auto ap = blade.g_a_prime();
        const auto omega = blade.parent_rotor.g_angular_velocity();
        const auto v_inf = blade.parent_rotor.g_wind_speed();
        const auto r = blade.g_radial_pos();
        const auto cs = blade.parent_rotor.speed_of_sound;

        const auto u = v_inf * ( 1 - a);
        const meters_per_second_t v = meters_per_second_t((omega * r * (1 + ap)).value());

        const auto w = math::hypot(u, v);
        const dimensionless_t mach = w / cs;

        if(mach <= 0.0 || mach >= 0.85) return; // outside valid area

        const auto beta = math::sqrt(math::max(1.0 - mach * mach, dimensionless_t(1e-12)));
        cl = cl / beta * (1 + (mach * mach) / (1 + beta));
    }



    inline void snel_stall_delay_correction(dimensionless_t& cl, dimensionless_t& cd, const BladeSection& blade, scalar_t K1 = 1.8, scalar_t K2 = 0.08) {
        using namespace units;
        using namespace units::math;

        const meter_t chord_len = blade.g_chord_len();
        const meter_t r = blade.g_radial_pos();

        const dimensionless_t fac = (chord_len / r) * (chord_len / r) * blade.g_local_tsr();
        const radian_t alpha = blade.g_alpha();
        if(alpha > 30.0_deg || alpha < -30.0_deg) return;

        const scalar_t s2 = math::sin(alpha) * math::sin(alpha);

        cl += K1 * fac * s2;
        cd += K2 * fac * s2;
    }


    inline void glauert_buhl_axial_induction(dimensionless_t& a, const dimensionless_t& F) {
    /*
     * From
     * "A New Empirical Relationship Between Thrust Coefficient and
     * Induction Factor for the Turbulent Windmill State” – NREL/TP-500-36834
    */
        if (a > 0.4) {
            const dimensionless_t CT = 4.0 * F * a * (1.0 - a);
            a = 0.5 * (2.0 + CT * (1.0 - 2.0/3.0) - math::sqrt((2.0 + CT * (1.0 - 2.0/3.0)) * (2.0 + CT * (1.0 - 2.0/3.0)) + 4.0 * (CT / 3.0 - 1.0)));
        }
    }

    inline void prandtl_tip_root(dimensionless_t& F, const BladeSection& blade, const radian_t& phi) {
        const dimensionless_t B = blade.parent_rotor.g_num_blades();
        const meter_t r = blade.g_radial_pos();
        const meter_t r_tip = blade.parent_rotor.g_rotor_radius();
        const meter_t r_hub = blade.parent_rotor.g_rotor_hub_radius();

        const dimensionless_t abs_phi = math::abs(math::sin(phi));

        if (abs_phi.value() < 1e-6) {
            F = dimensionless_t(1.0);
            return;
        }

        const dimensionless_t f_tip = 0.5 * B * ((r_tip - r) / r) / abs_phi;
        const dimensionless_t f_hub = 0.5 * B * ((r - r_hub) / r) / abs_phi;


        const dimensionless_t p = dimensionless_t(2.0 / M_PI);
        const dimensionless_t F_tip = p * dimensionless_t(math::acos(math::exp(-f_tip)).value());
        const dimensionless_t F_hub = p * dimensionless_t(math::acos(math::exp(-f_hub)).value());
        F = math::max(dimensionless_t(1e-4), F_tip * F_hub);
    }

}



