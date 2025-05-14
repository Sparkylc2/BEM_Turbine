#include "../include/headers.h"


// Example usage
int main() {
    // xf::demo();

    Rotor rotor("test");
    BladeSection test_section("2412", 0.5_m, 0.1_m, 0.1_rad, 0.01_m, rotor);
    rotor.run_bem();


    // test_section.update_alpha(radian_t(M_PI_4));
    // test_section.update_cl_cd();
    // std::cout << test_section.g_cl() << "\n";
    // std::cout << test_section.g_cd() << "\n";


    // rotor.save_rotor_json("test2_rotor");
    return 0;
}