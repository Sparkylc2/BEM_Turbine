//
// Created by Lukas Campbell on 10/05/2025.
//

#pragma once

#include "headers.h"
using namespace units;
using namespace units::literals;
using namespace units::length;
using namespace units::force;
using namespace units::angle;
using namespace units::dimensionless;
using namespace units::velocity;
using namespace units::mass;
using namespace units::volume;
using namespace units::density;
using namespace units::pressure;
using namespace units::time;

class BladeSection;

class Rotor {
    xf::Config cfg;

    std::string rotor_name;
    dimensionless_t num_blades = 0;
    meter_t rotor_radius = 0.0_m;
    meter_t rotor_hub_radius = 0.0_m;
    dimensionless_t tsr = 0.0;
    meters_per_second_t wind_speed = 0.0_mps;
    angular_velocity::radians_per_second_t angular_velocity = 0.0_rad_per_s;
    std::vector<BladeSection> blade_sections;


    newton_t drag = 0.0_N;
    newton_meter_t torque = 0.0_Nm;


    void initialize_rotor();
    nlohmann::json load_rotor_json();




public:
    using pascal_second_t = decltype(1.0_Pa * 1.0_s);

    xf::XFOILRunner runner;
    const kilograms_per_cubic_meter_t density = kilograms_per_cubic_meter_t(1.225);
    const pascal_second_t dynamic_viscosity = 1.813e-5_Pa * 1.0_s;
    const meters_per_second_t speed_of_sound = 343.0_mps;


    explicit Rotor(std::string name) : rotor_name{std::move(name)} { initialize_rotor(); }
    void save_rotor_json(const std::string& filename) const;
    void run_bem();




    const dimensionless_t& g_num_blades() const { return num_blades; }
    const meter_t& g_rotor_radius() const { return rotor_radius; }
    const dimensionless_t& g_rotor_tsr() const { return tsr; }
    const meters_per_second_t& g_wind_speed() const { return wind_speed; }
    const angular_velocity::radians_per_second_t& g_angular_velocity() const { return angular_velocity; }
    const meter_t& g_rotor_hub_radius() const { return rotor_hub_radius; }



};


