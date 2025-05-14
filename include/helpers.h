//
// Created by Lukas Campbell on 12/05/2025.
//

#pragma once
#include "headers.h"
namespace Helpers {

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
}
