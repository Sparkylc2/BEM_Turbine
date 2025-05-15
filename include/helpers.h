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

        const scalar_t sin_phi = math::copysign(math::max(math::fabs(cached_sin_phi), scalar_t(1e-3)), cached_sin_phi);
        const scalar_t cos_phi = math::copysign(math::max(math::fabs(cached_cos_phi), scalar_t(1e-3)), cached_cos_phi);

        return {sin_phi, cos_phi};
    }

    template <typename T>
    constexpr const T& clamp(const T& value, const T& lower, const T& upper) {
        if (value < lower) {
            return lower;
        } else if (value > upper) {
            return upper;
        }
        return value;
    }


    inline std::string sh_quote(const std::string& s) {
        std::string out{"'"};
        for (char c : s)
            if (c == '\'') out += "'\\''";   //   abc'd  ->  'abc'\''d'
            else out += c;
        out += '\'';
        return out;
    }

    inline std::string trim(const std::string& s) {
        const char* ws = " \t\r\n";
        const size_t start = s.find_first_not_of(ws);
        const size_t end   = s.find_last_not_of(ws);
        return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    }

    inline std::vector<std::string> split(const std::string& s) {
        std::istringstream iss(s);
        std::vector<std::string> tok;
        for(std::string w; iss >> w; ) tok.push_back(w);
        return tok;
    }

    inline std::string shell(const std::string& cmd, std::ostream* live = nullptr) {
        std::string out;
        char buf[256];

        #if defined(_WIN32)
                FILE* p = _popen(cmd.c_str(), "r");
        #else
                FILE* p = popen(cmd.c_str(), "r");
        #endif
                if(!p) return "";

                while(fgets(buf, sizeof(buf), p)) {
                    out += buf;
                    if(live) *live << buf;
                }

        #if defined(_WIN32)
                _pclose(p);
        #else
                pclose(p);
        #endif
                return out;
            }

            inline std::string get_project_root() {
        #ifdef XFOIL_PROJECT_ROOT
                return XFOIL_PROJECT_ROOT;
        #else
                throw std::runtime_error("XFOIL_PROJECT_ROOT macro not defined!");
        #endif
    }
}
