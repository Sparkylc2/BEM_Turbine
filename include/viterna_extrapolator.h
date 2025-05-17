// viterna_extrapolator.hpp  (unit-safe rewrite) ────────────────────────────
// Same algorithm but every literal and container now uses nholthaus/units
// types (`degree_t`, `dimensionless_t`) and units-aware math (`units::math`).
// No behaviour changed.  Include after <units/math.h> and your xf::Polar def.

#pragma once
#include "headers.h"        // brings in units + xf::Polar
#include <algorithm>
#include <vector>

using namespace units;
using namespace units::angle;
using namespace units::dimensionless;
using namespace units::literals;

class ViternaExtrapolator {
public:
    ViternaExtrapolator(const xf::Polar& polar, dimensionless_t cr75,
                        int n_alpha = 50) : cr75(cr75) {
        build_tables(polar, n_alpha);
    }

    std::pair<dimensionless_t, dimensionless_t> get_coefficients(degree_t alpha) const {
        alpha = math::fmod(alpha + 180.0_deg, 360.0_deg);
        if(alpha < 0.0_deg) alpha += 360.0_deg;
        alpha -= 180.0_deg;

        auto it = std::lower_bound(alphas.begin(), alphas.end(), alpha);
        if(it == alphas.begin()) return {c_ls.front(), c_ds.front()};
        if(it == alphas.end()) return {c_ls.back() , c_ds.back() };

        std::size_t i = std::distance(alphas.begin(), it);
        degree_t a0 = alphas[i - 1], a1 = alphas[i];
        auto t = ((alpha - a0) / (a1 - a0)).value();
        return {
            c_ls[i - 1] + (c_ls[i] - c_ls[i - 1]) * t,
            c_ds[i - 1] + (c_ds[i] - c_ds[i - 1]) * t
        };
    }

    bool is_in_original_range(degree_t a) const {
        return a >= alpha_min_orig && a <= alpha_max_orig;
    }

    const xf::Polar& get_extrapolated_polar() const { return extrapolated_polar; }

private:
    std::vector<degree_t> alphas;
    std::vector<dimensionless_t> c_ls, c_ds;
    xf::Polar extrapolated_polar;

    degree_t alpha_min_orig{0.0_deg}, alpha_max_orig{0.0_deg};
    dimensionless_t cr75;

    void build_tables(const xf::Polar& polar_in, int nalpha) {
        std::vector<degree_t> a_orig;
        std::vector<dimensionless_t> cl_orig, cd_orig;

        for(const auto& pt: polar_in.pts){
            a_orig .push_back(pt.alpha);
            cl_orig.push_back(pt.c_l);
            cd_orig.push_back(pt.c_d);
        }

        alpha_min_orig = a_orig.front();
        alpha_max_orig = a_orig.back();

        dimensionless_t AR = 1.0 / cr75;
        dimensionless_t cd_max  = std::max(*std::max_element(cd_orig.begin(), cd_orig.end()),
                                           1.11 + 0.018 * AR);

        auto it_ps    = std::max_element(cl_orig.begin(), cl_orig.end());
        int  i_ps = std::distance(cl_orig.begin(), it_ps);
        auto cl_ps = cl_orig[i_ps];
        auto cd_ps = cd_orig[i_ps];
        auto a_ps = a_orig[i_ps];

        int i_ns = 0;
        for(int i=0;i<i_ps;++i) if(cl_orig[i] < cl_orig[i_ns]) i_ns = i;
        auto cl_ns = cl_orig[i_ns];
        auto cd_ns = cd_orig[i_ns];
        auto a_ns = a_orig[i_ns];

        auto s=[](degree_t a){return math::sin(a);};
        auto c=[](degree_t a){return math::cos(a);};

        double B1pos = cd_max.value();
        double A1pos = 0.5*B1pos;
        double A2pos = (cl_ps - cd_max * s(a_ps) * c(a_ps)).value() * s(a_ps).value() / std::pow(c(a_ps).value(), 2);
        double B2pos = (cd_ps - cd_max*math::pow<2>(s(a_ps))).value() / c(a_ps).value();

        double B1neg = cd_max.value();
        double A1neg = 0.5 * B1neg;
        double A2neg = (cl_ns - cd_max*s(a_ns) * c(a_ns)).value() * s(a_ns).value() / std::pow(c(a_ns).value(),2);
        double B2neg = (cd_ns - cd_max*math::pow<2>(s(a_ns))).value() / c(a_ns).value();

        double da_pos = (alpha_max_orig - (-180.0_deg)).value() / nalpha;
        double da_neg = (180.0_deg + alpha_min_orig).value()  / nalpha;

        std::vector<degree_t> a_pos, a_neg;
        for(int i=0;i<nalpha;++i){
            a_pos.emplace_back(alpha_max_orig + degree_t((i+1)*da_pos));
            a_neg.emplace_back(alpha_min_orig - degree_t((i+1)*da_neg));
        }

        auto lift=[&](double A1,double A2,degree_t a){
            double s2a = math::sin(2.0*a).value();
            double s = math::sin(a).value();
            double c = math::cos(a).value();
            return dimensionless_t(A1 * s2a + A2 * c * c / std::max(1e-6,s));
        };

        auto drag=[&](double B1,double B2,degree_t a){
            double s = math::sin(a);
            double c = math::cos(a);
            return dimensionless_t(cd_max.value() * s * s + B2 * c);
        };

        std::vector<dimensionless_t> cl_pos,cd_pos,cl_neg,cd_neg;
        for(auto a: a_pos){ cl_pos.push_back(lift(A1pos,A2pos,a)); cd_pos.push_back(drag(B1pos,B2pos,a)); }
        for(auto a: a_neg){ cl_neg.push_back(lift(A1neg,A2neg,a)); cd_neg.push_back(drag(B1neg,B2neg,a)); }

        alphas.reserve(a_neg.size()+a_orig.size()+a_pos.size());
        c_ls.reserve(alphas.capacity()); c_ds.reserve(alphas.capacity());

        for(int i = static_cast<int>(a_neg.size()) - 1; i >= 0; --i) store_point(a_neg[i],cl_neg[i],cd_neg[i]);
        for(size_t i = 0; i < a_orig.size(); ++i) store_point(a_orig[i],cl_orig[i],cd_orig[i]);
        for(size_t i = 0; i < a_pos.size(); ++i) store_point(a_pos[i],cl_pos[i],cd_pos[i]);

        for(size_t i = 0; i < alphas.size(); ++i)
            extrapolated_polar.pts.push_back({alphas[i],c_ls[i],c_ds[i],dimensionless_t{0}});
    }

    void store_point(degree_t a, dimensionless_t cl, dimensionless_t cd){
        alphas.push_back(a); c_ls.push_back(cl); c_ds.push_back(cd);
    }
};