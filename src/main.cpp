#include "../include/headers.h"


using namespace matplot;
// Example usage
int main() {
    //
    // std::vector<double> time = Helpers::linspace(0, 3.0, 200);
    // std::vector<double> torque_vec;
    // std::vector<double> tsr_vec;
    //
    // double prev_tsr = 3.0;
    // double dt = time[1] - time[0];
    // for (int i = 0; i < time.size(); i++) {
    //     Rotor rotor = Rotor("blade_profile_test");
    //     rotor.initialize_rotor(dimensionless_t(prev_tsr));
    //
    //     double tsr = Rotor::simulate_rotor(rotor, dt);
    //
    //     std::cout << "Torque " << rotor.g_torque().value() << std::endl;
    //     torque_vec.push_back(rotor.g_torque().value());
    //     tsr_vec.push_back(tsr);
    //     std::cout << "TSR: " << tsr << std::endl;
    //     prev_tsr = tsr;
    // }
    //
    //
    // figure();
    // hold(on);
    // plot(time, tsr_vec) -> color("green");
    // plot(time, torque_vec) -> color("red");
    // xlabel("Time");
    // ylabel("TSR / Torque (Nm)");
    // title("TSR/Torque vs Time");
    // legend("TSR", "Torque");
    // grid(on);
    // show();







    double tsr_min = 0.05;
    double tsr_max = 8.0;

    std::vector<double> wind_speed_vec = Helpers::linspace(tsr_min, tsr_max, 200);
    std::vector<double> produced_power_vec;
    std::vector<double> cp_vec;
    std::vector<double> thrust_vec;
    std::vector<double> torque_vec;

    for (auto tsr : wind_speed_vec) {
        std::cout << "TSR: " << tsr << std::endl;
        Rotor rotor("blade_profile_test");
        rotor.initialize_rotor(dimensionless_t(tsr));
        rotor.run_bem();
        produced_power_vec.push_back(rotor.g_produced_power().value());
        cp_vec.push_back(rotor.g_c_p().value());
        thrust_vec.push_back(rotor.g_drag().value());
        torque_vec.push_back(rotor.g_torque().value());
    }

    double max_cp = *std::max_element(cp_vec.begin(), cp_vec.end());
    double max_cp_index = std::distance(cp_vec.begin(), std::max_element(cp_vec.begin(), cp_vec.end()));
    std::cout << "Max Cp: " << max_cp << " at TSR: " << wind_speed_vec[max_cp_index] << std::endl;
    std::cout << "Max Power: " << *std::max_element(produced_power_vec.begin(), produced_power_vec.end()) << std::endl;
    std::cout << "Max Thrust: " << *std::max_element(thrust_vec.begin(), thrust_vec.end()) << std::endl;
    std::cout << "Max Torque: " << *std::max_element(torque_vec.begin(), torque_vec.end()) << std::endl;


    // figure();
    // hold(on);
    // plot(wind_speed_vec, cp_vec) -> color("green");
    // xlabel("Tip-Speed Ratio (TSR)");
    // ylabel("C_p");
    // title("C_p vs TSR");
    // grid(on);
    // show();

    //
    //
    // figure();
    // hold(on);
    // auto l1 = plot(tsr_vec, produced_power_vec);
    // l1 -> color("blue");
    // l1 -> line_width(2.0);
    // l1 -> line_style("--");
    // l1 -> marker("o");
    // l1 -> marker_size(6.0);
    // plot(tsr_vec, total_power_vec)->color("red");
    // legend("Produced Power (, Total Power");
    // xlabel("Tip-Speed Ratio (TSR)");
    // ylabel("Power (W)");
    // title("Produced Power vs TSR");
    // grid(on);
    // show();
    //
    //
    // figure();
    // hold(on);
    // plot(tsr_vec, drag_vec)->color("blue");
    // plot(tsr_vec, torque_vec)->color("red");
    // legend("Drag (N), Torque Power (Nm)");
    // xlabel("Tip-Speed Ratio (TSR)");
    // ylabel("N/Nm");
    // title("Produced Power vs TSR");
    // grid(on);
    // show();
    //



    std::ostringstream file;
    file << (xf::Config::g_save_path("SAVE")) << "recent_run.csv";

    std::ofstream out(file.str());
    if (!out.is_open()) {
        std::cerr << "Could not open file " << file.str() << std::endl;
    }
    out << "TSR, THRUST (N), TORQUE (NM), PRODUCED POWER (W), C_P" << std::endl;

    for (size_t i = 0; i < wind_speed_vec.size(); i++) {
        out << wind_speed_vec[i] << ","
            << thrust_vec[i] << ","
            << torque_vec[i] << ","
            << produced_power_vec[i] << ","
            << cp_vec[i] << std::endl;
    }

    out.close();
    std::cout << "Saved results to  " << file.str() << std::endl;



    return 0;
}