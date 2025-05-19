// //
// // Created by Lukas Campbell on 12/05/2025.
// //
//
#pragma once
#include "headers.h"
namespace BEMSolver {


    void run_bem_solver(BladeSection& blade);
    std::pair<dimensionless_t, dimensionless_t> perform_iteration(const BladeSection& blade_ref, const dimensionless_t old_a, const dimensionless_t old_ap);
}

