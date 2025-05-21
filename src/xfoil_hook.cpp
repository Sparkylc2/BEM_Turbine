//
// Created by Lukas Campbell on 12/05/2025.
//
#include "../include/headers.h"

namespace xf {

    const std::unordered_map<std::string, std::string> Config::save_paths = {
        {"AIRFOIL", "/naca_data/airfoil_profiles/"},
        {"POLAR_TABLE", "/naca_data/polar_tables/"},
        {"BLADE", "/naca_data/blade_profiles/"},
        {"TMP", "/naca_data/run_tmp/"},
        {"SAVE", "/naca_data/bem_data/"}
    };
}