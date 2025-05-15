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
using namespace units::angular_velocity;
using namespace units::power;



void Rotor::initialize_rotor() {
    this -> runner = xf::XFOILRunner();
    this -> cfg = xf::Config();

    nlohmann::json j = load_rotor_json();

    this -> rotor_radius = meter_t(j["rotor_radius_m"].get<double>());
    this -> rotor_hub_radius = meter_t(j["rotor_hub_radius_m"].get<double>());
    this -> num_blades = dimensionless_t(j["num_blades"].get<int>());
    this -> tsr = dimensionless_t(j["tip_speed_ratio"].get<double>());
    this -> wind_speed = meters_per_second_t(j["wind_speed_mps"].get<double>());
    this -> angular_velocity = radians_per_second_t(((this -> tsr * this -> wind_speed) / this -> rotor_radius).value());

    this -> init_blade_sections(j);
}

void Rotor::initialize_rotor(double tsr) {
    this -> runner = xf::XFOILRunner();
    this -> cfg = xf::Config();

    nlohmann::json j = load_rotor_json();

    this -> rotor_radius = meter_t(j["rotor_radius_m"].get<double>());
    this -> rotor_hub_radius = meter_t(j["rotor_hub_radius_m"].get<double>());
    this -> num_blades = dimensionless_t(j["num_blades"].get<int>());
    this -> tsr = dimensionless_t(tsr);
    this -> wind_speed = meters_per_second_t(j["wind_speed_mps"].get<double>());
    this -> angular_velocity = radians_per_second_t(((this -> tsr * this -> wind_speed) / this -> rotor_radius).value());

    this -> init_blade_sections(j);
}


void Rotor::init_blade_sections(nlohmann::json& j) {

    for (size_t i = 0; i < j["blade"].size(); i++) {
        // double diff_r = j["blade"][i]["radial_pos_m"].get<double>() - j["blade"][i - 1]["radial_pos_m"].get<double>();
        meter_t diff_r = 0.0_m;
        if (i == 0) {
            diff_r = meter_t((j["blade"][i+1]["radial_pos_m"].get<double>() -
                      j["blade"][i]["radial_pos_m"].get<double>()) / 2.0 +
                      (j["blade"][i]["radial_pos_m"].get<double>() -
                       j["rotor_hub_radius_m"].get<double>()));
        } else if (i == j["blade"].size() - 1) {
            diff_r = meter_t((j["rotor_radius_m"].get<double>() -
                      j["blade"][i-1]["radial_pos_m"].get<double>()) / 2.0);
        } else {
            diff_r = meter_t((j["blade"][i+1]["radial_pos_m"].get<double>() -
                      j["blade"][i-1]["radial_pos_m"].get<double>()) / 2.0);
        }

        // if (i == 0 || i == j["blade"].size() - 1) diff_r /= 2.0;

        auto node = j["blade"][i];
        this -> blade_sections.emplace_back(
               node["naca"].get<std::string>(),
               meter_t(node["radial_pos_m"].get<double>()),
               meter_t(node["chord_len_m"].get<double>()),
               radian_t(node["twist_rad"].get<double>()),
               meter_t(diff_r),
               *this,
               node["coordinate_file"].get<std::string>()

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
    this -> drag = newton_t(0);
    this -> torque = newton_meter_t(0);

    for (auto& section : this -> blade_sections) {
        BEMSolver::run_bem_solver(section);
        section.post_bem_routine();

        if (!isnan(section.g_differential_drag().value())) this -> drag += section.g_differential_drag();
        if (!isnan(section.g_differential_torque().value())) this -> torque += section.g_differential_torque();
    }

    drag *= this -> num_blades;
    torque *= this -> num_blades;

    produced_power = watt_t((this -> torque * this -> angular_velocity).value());
    total_power = 0.5 * M_PI * density * rotor_radius * rotor_radius * wind_speed * wind_speed * wind_speed;

    cp = produced_power / total_power;
    //
    // std::cout << "Calculated Axial Drag: " << this->drag << "\n";
    // std::cout << "Calculated Torque: " << this->torque << "\n";
    // std::cout << "Mechanical Power: " << produced_power << " W\n";
    // std::cout << "Available Wind Power: " << total_power << " W\n";
    // std::cout << "Power Coefficient (Cp): " << cp << " [dimensionless]\n";

    if (cp.value() > 0.593) {
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
            {"radial_pos_m", section.g_radial_pos().value()},
            {"chord_len_m",  section.g_chord_len().value()},
            {"twist_rad",    section.g_twist_angle().value()}
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

