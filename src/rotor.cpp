//
// Created by Lukas Campbell on 10/05/2025.
//

#include "../include/headers.h"

using namespace units;
using namespace units::literals;
using namespace units::length;
using namespace units::force;
using namespace units::angle;
using namespace units::dimensionless;
using namespace units::velocity;


void Rotor::initialize_rotor() {
    this -> runner = xf::XFOILRunner();
    this -> cfg = xf::Config();

    nlohmann::json j = load_rotor_json();

    this -> rotor_radius = j["rotor_radius_m"].get<double>() * meter_t(1.0);
    this -> rotor_hub_radius = j["rotor_hub_radius_m"].get<double>() * meter_t(1.0);
    this -> num_blades = j["num_blades"].get<int>() * dimensionless_t(1.0);
    this -> tsr = j["tip_speed_ratio"].get<double>() * dimensionless_t(1.0);
    this -> wind_speed = j["wind_speed_mps"].get<double>() * meters_per_second_t(1.0);
    this -> angular_velocity = angular_velocity::radians_per_second_t(((this -> tsr * this -> wind_speed) / this -> rotor_radius).value());

    for (size_t i = 0; i < j["blade"].size(); i++) {
        // double diff_r = j["blade"][i]["radial_pos_m"].get<double>() - j["blade"][i - 1]["radial_pos_m"].get<double>();
        double diff_r = 0.0;
        if (i == 0) {
            diff_r = (j["blade"][i+1]["radial_pos_m"].get<double>() -
                      j["blade"][i]["radial_pos_m"].get<double>()) / 2.0 +
                      (j["blade"][i]["radial_pos_m"].get<double>() -
                       j["rotor_hub_radius_m"].get<double>());
        } else if (i == j["blade"].size() - 1) {
            diff_r = (j["rotor_radius_m"].get<double>() -
                      j["blade"][i-1]["radial_pos_m"].get<double>()) / 2.0;
        } else {
            diff_r = (j["blade"][i+1]["radial_pos_m"].get<double>() -
                      j["blade"][i-1]["radial_pos_m"].get<double>()) / 2.0;
        }

        // if (i == 0 || i == j["blade"].size() - 1) diff_r /= 2.0;

        auto node = j["blade"][i];
        this -> blade_sections.emplace_back(
               node["naca"].get<std::string>(),
               node["radial_pos_m"].get<double>() * meter_t(1.0),
               node["chord_len_m"].get<double>() * meter_t(1.0),
               node["twist_rad"].get<double>() * radian_t(1.0),
               meter_t(diff_r),
               *this
       );
    }
}



nlohmann::json Rotor::load_rotor_json() {
    std::ifstream rotor_file;

    if (rotor_name.empty()) std::cerr << "Rotor name is empty. Cannot load rotor data.\n";

    std::string name = xf::Config::g_save_path("BLADE") + rotor_name + ".json";

    rotor_file.open(name);

    if (!rotor_file.is_open()) std::cerr << "Failed to open rotor file: " << name << "\n";

    nlohmann::json j;
    rotor_file >> j;
    return j;
}

void Rotor::run_bem() {
    this->drag = newton_t(0);
    this->torque = newton_meter_t(0);

    for (auto& section : this -> blade_sections) {
        BEMSolver::run_bem_solver(section, true);
        section.post_bem_routine();
        this -> drag += section.g_differential_drag();
        this -> torque += section.g_differential_torque();
    }

    drag *= this -> num_blades;
    torque *= this -> num_blades;

    auto produced_power = this -> torque * this -> angular_velocity;
    auto total_power = 0.5 * M_PI * density * rotor_radius * rotor_radius * wind_speed * wind_speed * wind_speed;

    auto power_coefficient = produced_power / total_power;

    std::cout << "Calculated Axial Drag: " << this->drag << "\n";
    std::cout << "Calculated Torque: " << this->torque << "\n";
    std::cout << "Mechanical Power: " << produced_power << " W\n";
    std::cout << "Available Wind Power: " << total_power << " W\n";
    std::cout << "Power Coefficient (Cp): " << power_coefficient << " [dimensionless]\n";

    if (power_coefficient.value() > 0.593) {
        std::cout << "WARNING: Calculated cp exceeds betz limit\n";
    }
}




void Rotor::save_rotor_json(const std::string& name) const {
    nlohmann::json j;

    j["rotor_hub_radius_m"] = unit_cast<double>(rotor_hub_radius);
    j["rotor_radius_m"]  = unit_cast<double>(rotor_radius);
    j["num_blades"]      = unit_cast<int>(num_blades);
    j["tip_speed_ratio"] = unit_cast<double>(tsr);
    j["wind_speed_mps"]  = unit_cast<double>(wind_speed);

    j["blade"] = nlohmann::json::array();

    for (const auto& section : this->blade_sections) {
        j["blade"].push_back({
            {"naca",         section.g_naca_code()},
            {"radial_pos_m", unit_cast<double>(section.g_radial_pos())},
            {"chord_len_m",  unit_cast<double>(section.g_chord_len())},
            {"twist_rad",    unit_cast<double>(section.g_twist_angle())}
        });
    }

    const std::string path = xf::Config::g_save_path("BLADE", name + ".json");

    try {
        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
        std::ofstream out(path);
        if (!out.is_open()) {
            std::cerr << "Failed to open output file: " << path << "\n";
            return;
        }

        out << std::setw(4) << j << "\n";
        std::cout << "Saved rotor configuration to: " << path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error writing rotor JSON: " << e.what() << "\n";
    }
}

