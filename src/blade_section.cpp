//
// Created by Lukas Campbell on 10/05/2025.
//

#include "../include/headers.h"

void BladeSection::initialize_blade_section() {
    local_solidity = parent_rotor.g_num_blades() * chord_length / (2.0 * M_PI * radial_position);
    local_tip_speed_ratio = meters_per_second_t((parent_rotor.g_angular_velocity() * radial_position).value()) / parent_rotor.g_wind_speed();

    non_dim_radius = radial_position / parent_rotor.g_rotor_radius();

    angular_velocity = parent_rotor.g_angular_velocity();

    double W = std::sqrt(parent_rotor.g_wind_speed().value() * parent_rotor.g_wind_speed().value() + (angular_velocity * radial_position).value() * (angular_velocity * radial_position).value());
    double re = parent_rotor.density.value() * W * chord_length.value() / parent_rotor.dynamic_viscosity.value();
    double mach = W / parent_rotor.speed_of_sound.value();



    cfg = xf::Config();
    cfg.naca = naca_code;
    if (!coordinate_file.empty()) {
        if (parent_rotor.runner.get_runtime() == "docker") {
            cfg.coord_file = "/work/naca_data/airfoil_profiles/" + coordinate_file;
            // std::cout << "Using coordinate file: " << cfg.coord_file << "\n";
        } else {
            cfg.coord_file = "../naca_data/airfoil_profiles/" + coordinate_file;
        }
    }
    cfg.alpha = 0.0_deg;
    cfg.re = re;
    cfg.mach = mach;




    a = dimensionless_t(1.0 / 3.0);
    a_prime = dimensionless_t(1.0 / 3.0);
}





void BladeSection::update_alpha(radian_t phi) {
    this -> alpha = phi - twist_angle;
    cfg.alpha = alpha;
}


void BladeSection::update_cl_cd() {
    static std::unordered_map<std::string, xf::Polar> basePolars;

    std::string re = std::to_string(static_cast<int>(cfg.re.value()));
    xf::Polar basePolar;
    std::string polarName = naca_code.empty() ? coordinate_file + re: naca_code + re;
    auto it = basePolars.find(polarName);

    if (it == basePolars.end()) {
        std::cout << "Calculating polar for: " << polarName << std::endl;
        xf::XFOILRunner &runner = parent_rotor.runner;

        auto routine = [&](std::ofstream &ss) {
            using namespace xf;
            XFOILRunner::toggle_headless(ss);
            XFOILRunner::init_airfoil_geometry(ss, cfg);
            XFOILRunner::configure_airfoil(ss, cfg);
            XFOILRunner::configure_solver(ss, cfg);
            XFOILRunner::set_alpha_sweep(ss, cfg);
            XFOILRunner::save_airfoil(ss, cfg);
            XFOILRunner::quit(ss);
        };

        try {
            if (!runner.run(cfg, basePolar, routine) || basePolar.pts.empty()) {
                if (naca_code.empty()) {
                    // std::cerr << "[XFOIL] failed for airfoil file: " << coordinate_file << "\n";
                // } else {
                    // std::cerr << "[XFOIL] failed for airfoil: " << naca_code << "\n";
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "XFOIL exception: " << e.what() << std::endl;
        }

        basePolars.emplace(polarName, basePolar);
    } else {
        basePolar = it -> second;
    }

    const dimensionless_t cr75 = chord_length / (parent_rotor.g_rotor_radius() * 0.75);

    ViternaExtrapolator extrap(basePolar, cr75);

    auto [cl_val, cd_val]  = extrap.get_coefficients(alpha);


    if (std::isnan(cl_val.value()) || std::isnan(cd_val.value())) {
        std::cerr << "XFOIL failed to generate points for airfoil: " + (naca_code.empty() ? cfg.coord_file : naca_code) << std::endl;
    }

    // if (!extrap.is_in_original_range(alpha))

    c_l = cl_val;
    c_d = cd_val;
}


void BladeSection::update_induction_factors(const dimensionless_t a, const dimensionless_t a_prime) {
    this -> a = a;
    this -> a_prime = a_prime;
}

void BladeSection::post_bem_routine() {
    // std::cout << "Running post BEM routine for blade section " << naca_code << " at radius: " << radial_position << "\n";
    const radian_t inflow_angle = alpha + twist_angle;

    meters_per_second_t w1 = parent_rotor.g_wind_speed() * (1.0 - a) / math::sin(inflow_angle);
    meters_per_second_t w2 = meters_per_second_t((parent_rotor.g_angular_velocity() * radial_position).value()) * (1.0 + a_prime) / math::cos(inflow_angle);

    if (math::abs(w1 - w2) > 1e-6_mps) {
        std::cerr << "Warning: Inflow angle is not converged!  " << "Residual = " << w1 - w2 << "\n";
        w1 = 0_mps;
        w2 = 0_mps;
    }
    this -> w = (w1 + w2) / 2.0;
    compute_differential_forces();
}


void BladeSection::compute_differential_forces() {
    const radian_t phi = alpha + twist_angle;
    update_cl_cd();
    const dimensionless_t c_n = c_l * math::cos(phi) + c_d * math::sin(phi);
    const dimensionless_t c_t = c_l * math::sin(phi) - c_d * math::cos(phi);


    differential_drag   = 0.5 * parent_rotor.density * w * w * c_n * chord_length * differential_radius;
    differential_torque = 0.5 * parent_rotor.density * w * w * c_t * chord_length * radial_position * differential_radius;
}


void BladeSection::update_reynolds(const meters_per_second_t& w) {
    cfg.re = parent_rotor.density * w * chord_length / parent_rotor.dynamic_viscosity;
}

void BladeSection::update_mach(const meters_per_second_t &w) {
    cfg.mach = w / parent_rotor.speed_of_sound;
}



