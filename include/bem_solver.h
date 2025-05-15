//
// Created by Lukas Campbell on 12/05/2025.
//

#pragma once
#include "headers.h"
namespace BEMSolver {
    bool bracket(const std::function<double(BladeSection&, double x)> &f, BladeSection& blade, double min, double max, int n, bool backwards_search, double& b_min, double& b_max);

    void run_bem_solver(BladeSection& blade);
    void perform_blade_iteration(const BladeSection& blade, radian_t& phi);

    dimensionless_t compute_residual(const BladeSection& blade_ref, const radian_t& phi);

    double newton_raphson(
        const std::function<double(BladeSection&, double)>&  f,
        const std::function<double(BladeSection&, double)>&  df,
        BladeSection&                                        blade,
        double                                               x0,
        double                                               tol,
        int                                                  max_iter);

    double numerical_derivative(const BladeSection& blade,
                                double              x,
                                double              h);

    double bisection_method(const std::function<double(BladeSection&, double)> &f,
                           BladeSection& blade, double a, double b,
                           const double tol, const int max_iter);




}

