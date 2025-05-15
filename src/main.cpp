#include "../include/headers.h"


using namespace matplot;
// Example usage
int main() {
    // xf::demo();

    double tsr_min = 1.0;
    double tsr_max = 4.5;

    std::vector<double> tsr_vec = Helpers::linspace(tsr_min, tsr_max, 50);

    std::vector<double> produced_power_vec, total_power_vec;
    std::vector<double> cp_vec;
    std::vector<double> torque_vec;
    std::vector<double> drag_vec;

    for (auto tsr: tsr_vec) {
        std::cout << "TSR: " << tsr << std::endl;
        Rotor rotor("blade a");
        rotor.initialize_rotor(tsr);
        rotor.run_bem();
        produced_power_vec.push_back(rotor.g_produced_power().value());
        total_power_vec.push_back(rotor.g_total_power().value());
        cp_vec.push_back(rotor.g_c_p().value());
        torque_vec.push_back(rotor.g_torque().value());
        drag_vec.push_back(rotor.g_drag().value());
    }




    figure();
    hold(on);

    // plot(tsr_vec, produced_power_vec)->color("blue");
    // plot(tsr_vec, total_power_vec)->color("red");
    plot(tsr_vec, cp_vec)->color("green");

    legend("Produced Power");
    xlabel("Tip-Speed Ratio (TSR)");
    ylabel("Power / Cp");
    title("BEM Results vs TSR");
    grid(on);

    show();

    return 0;
}