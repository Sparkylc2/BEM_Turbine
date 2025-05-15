//
// Created by Lukas Campbell on 12/05/2025.
//

#pragma once
#include "headers.h"
namespace Helpers {
    using namespace units;
    using namespace units::angle;
    using namespace units::dimensionless;

    inline std::vector<double> linspace(const double& min, const double& max, const int& num_p) {
        std::vector<double> vec;
        if (num_p <= 0) std::cerr << "Error: num_p must be greater than 0\n" << std::endl;

        if (num_p == 1) {
            vec.push_back(min);
            return vec;
        }

        const double step = (max - min) / (num_p - 1);
        for (int i = 0; i < num_p; ++i) {
            vec.push_back(min + i * step);
        }
        return vec;

    }

    inline bool is_approx(double a, double b, double eps) {
        return std::abs(a - b) < eps;
    }


    inline std::pair<scalar_t, scalar_t> get_clamped_phi_components(const radian_t& phi) {
        const scalar_t cached_sin_phi = math::sin(phi);
        const scalar_t cached_cos_phi = math::cos(phi);

        const scalar_t sin_phi = math::copysign(math::max(math::fabs(cached_sin_phi), 1e-3), cached_sin_phi);
        const scalar_t cos_phi = math::copysign(math::max(math::fabs(cached_cos_phi), 1e-3), cached_cos_phi);

        return {sin_phi, cos_phi};
    }

}
