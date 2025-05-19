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




    double tsr_min = 0.0;
    double tsr_max = 12;

    std::vector<double> tsr_vec = Helpers::linspace(tsr_min, tsr_max, 200);

    std::vector<double> produced_power_vec, total_power_vec;
    std::vector<double> cp_vec;
    std::vector<double> torque_vec;
    std::vector<double> drag_vec;

    for (auto tsr: tsr_vec) {
        std::cout << "TSR: " << tsr << std::endl;
        Rotor rotor("NREL_5MW_Blade");
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
    plot(tsr_vec, cp_vec) -> color("green");
    xlabel("Tip-Speed Ratio (TSR)");
    ylabel("C_p");
    title("C_p vs TSR");
    grid(on);
    show();


    figure();
    hold(on);
    auto l1 = plot(tsr_vec, produced_power_vec);
    l1 -> color("blue");
    l1 -> line_width(2.0);
    l1 -> line_style("--");
    l1 -> marker("o");
    l1 -> marker_size(6.0);
    plot(tsr_vec, total_power_vec)->color("red");
    legend("Produced Power (, Total Power");
    xlabel("Tip-Speed Ratio (TSR)");
    ylabel("Power (W)");
    title("Produced Power vs TSR");
    grid(on);
    show();


    figure();
    hold(on);
    plot(tsr_vec, drag_vec)->color("blue");
    plot(tsr_vec, torque_vec)->color("red");
    legend("Drag (N), Torque Power (Nm)");
    xlabel("Tip-Speed Ratio (TSR)");
    ylabel("N/Nm");
    title("Produced Power vs TSR");
    grid(on);
    show();


    return 0;
}