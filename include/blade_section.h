//
// Created by Lukas Campbell on 10/05/2025.
//

#pragma once

using namespace units;
using namespace units::literals;
using namespace units::length;
using namespace units::force;
using namespace units::angle;
using namespace units::dimensionless;
using namespace units::torque;
using namespace units::velocity;

class Rotor;

class BladeSection {

    xf::Config cfg;

    std::string naca_code;
    std::string coordinate_file;


    meter_t radial_position;
    dimensionless_t non_dim_radius;
    radian_t angular_position;
    angular_velocity::radians_per_second_t angular_velocity;

    meter_t chord_length;
    radian_t twist_angle;
    dimensionless_t local_tip_speed_ratio;
    dimensionless_t local_solidity;

    dimensionless_t c_l;
    dimensionless_t c_d;

    dimensionless_t c_n;
    dimensionless_t c_t;

    radian_t alpha;

    dimensionless_t a;
    dimensionless_t a_prime;

    meters_per_second_t w;

    newton_t differential_drag;
    newton_meter_t differential_torque;

    meter_t differential_radius;

    void initialize_blade_section();
public:
    explicit BladeSection(const std::string &naca_code,
        const meter_t radial_position,
        const meter_t chord_length,
        const radian_t twist_angle,
        const meter_t differential_radius,
        Rotor& parent_rotor,
        const std::string coordinate_file = ""
    )
    :
        naca_code{naca_code},
        radial_position{radial_position},
        chord_length{chord_length},
        twist_angle{twist_angle},
        differential_radius{differential_radius},
        coordinate_file{coordinate_file},
        parent_rotor{parent_rotor}

    {
        initialize_blade_section();
    }



    Rotor& parent_rotor;

    void compute_differential_forces();

    void update_alpha(radian_t alpha);
    void update_reynolds(const meters_per_second_t& w);
    void update_mach(const meters_per_second_t& w);

    void update_cl_cd();
    void post_bem_routine();

    void update_induction_factors(dimensionless_t a, dimensionless_t a_prime);



    void set_differential_drag(newton_t dD) { this -> differential_drag = dD; }
    void set_differential_torque(newton_meter_t dT) { this -> differential_torque = dT; }

    // getters
    [[nodiscard]] meter_t const& g_chord_len() const { return chord_length; }
    [[nodiscard]] meter_t const& g_radial_pos() const { return radial_position; }
    [[nodiscard]] radian_t const& g_angular_pos() const { return angular_position; }
    [[nodiscard]] angular_velocity::radians_per_second_t const& g_angular_vel() const { return angular_velocity; }
    [[nodiscard]] dimensionless_t const& g_non_dim_radius() const { return non_dim_radius; }
    [[nodiscard]] radian_t const& g_twist_angle() const { return twist_angle; }
    [[nodiscard]] dimensionless_t const& g_local_tsr() const { return local_tip_speed_ratio; }
    [[nodiscard]] std::string const& g_naca_code() const { return naca_code; }
    [[nodiscard]] dimensionless_t const& g_local_solidity() const { return local_solidity; }
    [[nodiscard]] dimensionless_t const& g_cl() const { return c_l; }
    [[nodiscard]] dimensionless_t const& g_cd() const { return c_d; }
    [[nodiscard]] radian_t const& g_alpha() const { return alpha; }
    [[nodiscard]] dimensionless_t const& g_a() const { return a; }
    [[nodiscard]] dimensionless_t const& g_a_prime() const { return a_prime; }
    [[nodiscard]] dimensionless_t const& g_c_n() const { return c_n; }
    [[nodiscard]] dimensionless_t const& g_c_t() const { return c_t; }
    [[nodiscard]] meters_per_second_t const& g_w() const { return w; }
    [[nodiscard]] newton_t const& g_differential_drag() const { return differential_drag; }
    [[nodiscard]] newton_meter_t const& g_differential_torque() const { return differential_torque; }
    [[nodiscard]] meter_t const& g_differential_radius() const { return differential_radius; }


};



