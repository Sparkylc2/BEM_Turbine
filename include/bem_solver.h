//
// Created by Lukas Campbell on 12/05/2025.
//

#pragma once
#include "headers.h"
namespace BEMSolver {
    using namespace units;
    using namespace units::literals;
    using namespace units::length;
    using namespace units::force;
    using namespace units::angle;
    using namespace units::dimensionless;



    dimensionless_t tip_hub_loss_correction(meter_t r, meter_t r_hub, meter_t r_tip, radian_t inflow_angle, dimensionless_t num_blades);
    void glauert_prandtl_mach_correction(dimensionless_t& cl, dimensionless_t& cd, dimensionless_t& mach);
    void reynolds_number_correction(dimensionless_t& cl, dimensionless_t& cd, dimensionless_t& re);



    bool bracket(const std::function<double(BladeSection&, double x)> &f, BladeSection& blade, double min, double max, int n, bool backwards_search, double& b_min, double& b_max);

    void run_bem_solver(BladeSection& blade, bool eps_all);
    dimensionless_t compute_residual(const BladeSection& blade, const radian_t& phi);

    double newton_raphson(
        const std::function<double(BladeSection&, double)>&  f,
        const std::function<double(BladeSection&, double)>&  df,
        BladeSection&                                        blade,
        double                                               x0,
        const double                                         tol,
        const int                                            max_iter);

    double numerical_derivative(const BladeSection& blade, double x, double h);

    double bisection_method(const std::function<double(BladeSection&, double)> &f,
                           BladeSection& blade, double a, double b,
                           const double tol, const int max_iter);
    double brent(const std::function<double(BladeSection&, double)> &f,
                      BladeSection& blade, double a, double b,
                      const double tol, const int max_iter);




}

