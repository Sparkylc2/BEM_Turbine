//
// Created by Lukas Campbell on 12/05/2025.
//

#pragma once
#include "headers.h"
namespace BEMSolver {
    using Blade_Func = std::function<dimensionless_t(BladeSection&, const radian_t&)>;
    radian_t newton_raphson(const Blade_Func& f, const Blade_Func& df, BladeSection& blade, radian_t phi_0, scalar_t tol, int max_iter);
    double numerical_derivative(const BladeSection& blade, double x, double h);


    void run_bem_solver(BladeSection& blade);
    std::pair<dimensionless_t, dimensionless_t> perform_blade_iteration(BladeSection &blade, const radian_t& phi);
    radian_t solve_phi_with_retries(BladeSection& blade, Blade_Func& f, Blade_Func& df);
    dimensionless_t compute_residual(const BladeSection& blade_ref, const radian_t& phi);
}

