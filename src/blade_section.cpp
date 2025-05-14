//
// Created by Lukas Campbell on 10/05/2025.
//

#include "../include/headers.h"

void BladeSection::initialize_blade_section() {
    local_solidity = parent_rotor.g_num_blades() * chord_length / (2.0 * M_PI * radial_position);
    local_tip_speed_ratio = dimensionless_t((parent_rotor.g_angular_velocity() * radial_position / parent_rotor.g_wind_speed()).value());

    non_dim_radius = radial_position / parent_rotor.g_rotor_radius();

    angular_velocity = parent_rotor.g_angular_velocity();

    cfg = xf::Config();
    cfg.naca = naca_code;
    cfg.alpha = 0.0;
}




void BladeSection::update_alpha(radian_t phi) {
    this -> alpha = phi - twist_angle;
    cfg.alpha = degree_t(alpha).value();
}


void BladeSection::update_cl_cd()
{
    static std::unordered_map<std::string, xf::Polar> basePolars;

    xf::Polar basePolar;
    auto it = basePolars.find(naca_code);
    if (it == basePolars.end())
    {
        xf::XFOILRunner &runner = parent_rotor.runner;

        auto routine = [&](std::ofstream &ss)
        {
            using namespace xf;
            XFOILRunner::toggle_headless(ss);
            XFOILRunner::init_airfoil_geometry(ss, cfg);
            XFOILRunner::configure_airfoil(ss, cfg);
            XFOILRunner::configure_solver(ss, cfg);
            XFOILRunner::set_alpha_sweep(ss, cfg);
            XFOILRunner::quit(ss);
        };


        if (!runner.run(cfg, basePolar, routine) || basePolar.pts.empty())
        {
            std::cerr << "[XFOIL]   failed for airfoil " << naca_code
                      << " – using flat-plate model\n";
        }

        basePolars.emplace(naca_code, basePolar);
    }
    else
        basePolar = it->second;

    const double cr75 = chord_length.value() /
                        (parent_rotor.g_rotor_radius().value() * 0.75);

    ViternaExtrapolator extrap(basePolar, cr75);

    const double alpha_deg = degree_t(alpha).value();
    auto [cl_val, cd_val]  = extrap.getCoefficients(alpha_deg);


    if (std::isnan(cl_val) || std::isnan(cd_val)) {
        throw std::runtime_error("XFOIL failed to generate points for Re=" + std::to_string(cfg.re));
    }

    if (!extrap.isInOriginalRange(alpha_deg))

    cl = dimensionless_t(cl_val);
    cd = dimensionless_t(cd_val);
}


void BladeSection::update_induction_factors(dimensionless_t a, dimensionless_t a_prime) {
    this -> a = a;
    this -> a_prime = a_prime;
}

void BladeSection::post_bem_routine() {
    std::cout << "Running post BEM routine for blade section " << naca_code << " at radius: " << unit_cast<double>(radial_position) << " m\n";
    const radian_t inflow_angle = alpha + twist_angle;

    const double w1 = unit_cast<double>(parent_rotor.g_wind_speed() * (1.0 - a.value()) / math::sin(inflow_angle));
    const double w2 = unit_cast<double>(parent_rotor.g_angular_velocity() * radial_position * (1.0 + a_prime.value()) / math::cos(inflow_angle));
    if (abs(w1 - w2) > 0.01) {
        std::cerr << "Warning: Inflow angle is not converged!\n";
    }
    this -> w = meters_per_second_t(w1);

    compute_differential_forces();
}


void BladeSection::compute_differential_forces() {
    const radian_t inflow_angle = alpha + twist_angle;
    const dimensionless_t c_n = cl * math::cos(inflow_angle) - cd * math::sin(inflow_angle);
    const dimensionless_t c_t = cl * math::sin(inflow_angle) + cd * math::cos(inflow_angle);

    differential_drag = 0.5 * parent_rotor.density * w * w * c_n * chord_length * differential_radius;
    differential_torque = 0.5 * parent_rotor.density * w * w * c_t * chord_length * radial_position * differential_radius;
}


void BladeSection::update_reynolds(const meters_per_second_t& w) {
    cfg.re = parent_rotor.density * w * chord_length / parent_rotor.dynamic_viscosity;
}

void BladeSection::update_mach(const meters_per_second_t &w) {
    cfg.mach = w / parent_rotor.speed_of_sound;
}



