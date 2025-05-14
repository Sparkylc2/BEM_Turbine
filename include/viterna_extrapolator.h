// ============================================================================
//  ViternaExtrapolator.hpp / .cpp   ―   fully validated implementation
//  Formulation: Viterna & Corrigan 1981, NASA TM-81182
// ============================================================================

class ViternaExtrapolator
{
public:
    ViternaExtrapolator(const xf::Polar &polar, double cr75, int nalpha = 50)
        : cr75_(cr75)
    {
        buildTables(polar, nalpha);
    }

    std::pair<double,double> getCoefficients(double alpha_deg) const
    {
        alpha_deg = std::fmod(alpha_deg + 180.0, 360.0);
        if (alpha_deg < 0.0) alpha_deg += 360.0;
        alpha_deg -= 180.0;

        const auto it = std::lower_bound(alphas_.begin(), alphas_.end(), alpha_deg);
        if (it == alphas_.begin())          return { cls_.front(), cds_.front() };
        if (it == alphas_.end())            return { cls_.back(),  cds_.back()  };

        std::size_t i = std::distance(alphas_.begin(), it);
        const double a0 = alphas_[i-1], a1 = alphas_[i];
        const double t  = (alpha_deg - a0) / (a1 - a0);

        return { cls_[i-1] + t * (cls_[i] - cls_[i-1]),
                 cds_[i-1] + t * (cds_[i] - cds_[i-1]) };
    }

    bool isInOriginalRange(double a_deg) const
    {
        return  a_deg >= alpha_min_orig_  &&  a_deg <= alpha_max_orig_;
    }

    const xf::Polar &getExtrapolatedPolar() const { return polar_ex_; }

private:
    std::vector<double> alphas_, cls_, cds_;
    xf::Polar           polar_ex_;
    double              alpha_min_orig_ = 0.0;
    double              alpha_max_orig_ = 0.0;
    double              cr75_;

    static constexpr double deg2rad = M_PI / 180.0;


    void buildTables(const xf::Polar &polar_in, int nalpha)
    {

        std::vector<double> a_orig, cl_orig, cd_orig;
        for (const auto &pt : polar_in.pts) {
            a_orig.push_back(pt.alpha);
            cl_orig.push_back(pt.cl);
            cd_orig.push_back(pt.cd);
        }
        alpha_min_orig_ = a_orig.front();
        alpha_max_orig_ = a_orig.back();


        const double AR      = 1.0 / cr75_;
        const double cd_max  = std::max(*std::max_element(cd_orig.begin(), cd_orig.end()),
                                        1.11 + 0.018 * AR);


        const auto it_ps     = std::max_element(cl_orig.begin(), cl_orig.end());
        const int   i_ps     = std::distance(cl_orig.begin(), it_ps);
        const double cl_ps   = cl_orig[i_ps];
        const double cd_ps   = cd_orig[i_ps];
        const double a_ps    = a_orig [i_ps];


        int i_ns = 0;
        for (int i = 0; i < i_ps; ++i)
            if (cl_orig[i] < cl_orig[i_ns]) i_ns = i;

        const double cl_ns = cl_orig[i_ns];
        const double cd_ns = cd_orig[i_ns];
        const double a_ns  = a_orig [i_ns];


        auto s = [](double a_deg){ return std::sin(a_deg * deg2rad); };
        auto c = [](double a_deg){ return std::cos(a_deg * deg2rad); };


        const double B1pos = cd_max;
        const double A1pos = 0.5*B1pos;
        const double A2pos = (cl_ps - cd_max * s(a_ps)*c(a_ps))
                             * s(a_ps) / (c(a_ps)*c(a_ps));
        const double B2pos = (cd_ps - cd_max * s(a_ps)*s(a_ps)) / c(a_ps);

        const double B1neg = cd_max;
        const double A1neg = 0.5*B1neg;
        const double A2neg = (cl_ns - cd_max * s(a_ns)*c(a_ns))
                             * s(a_ns) / (c(a_ns)*c(a_ns));
        const double B2neg = (cd_ns - cd_max * s(a_ns)*s(a_ns)) / c(a_ns);


        const double da_pos = ( 180.0 - alpha_max_orig_) / nalpha;
        const double da_neg = ( alpha_min_orig_ + 180.0) / nalpha;

        std::vector<double> a_pos, a_neg;
        for (int i = 0; i < nalpha; ++i) {
            a_pos.push_back(alpha_max_orig_ + (i+1)*da_pos);
            a_neg.push_back(alpha_min_orig_ - (i+1)*da_neg);
        }

        auto lift = [](double A1, double A2, double a_deg){
            const double s2a = std::sin(2.0 * a_deg * deg2rad);
            const double  s  = std::sin(    a_deg * deg2rad);
            const double  c  = std::cos(    a_deg * deg2rad);
            return A1 * s2a + A2 * (c*c) / std::max(1e-6, s);
        };
        auto drag = [cd_max](double B1, double B2, double a_deg){
            const double s = std::sin(a_deg * deg2rad);
            const double c = std::cos(a_deg * deg2rad);
            return cd_max * s*s + B2 * c;
        };

        std::vector<double> cl_pos, cd_pos, cl_neg, cd_neg;
        for (double a : a_pos) {
            cl_pos.push_back(lift(A1pos, A2pos, a));
            cd_pos.push_back(drag(B1pos, B2pos, a));
        }
        for (double a : a_neg) {
            cl_neg.push_back(lift(A1neg, A2neg, a));
            cd_neg.push_back(drag(B1neg, B2neg, a));
        }


        alphas_.reserve(a_neg.size() + a_orig.size() + a_pos.size());
        cls_.reserve(alphas_.capacity());
        cds_.reserve(alphas_.capacity());


        for (int i = (int)a_neg.size()-1; i >= 0; --i){
            storePoint(a_neg[i], cl_neg[i], cd_neg[i]);
        }
        for (std::size_t i = 0; i < a_orig.size(); ++i){
            storePoint(a_orig[i], cl_orig[i], cd_orig[i]);
        }
        for (std::size_t i = 0; i < a_pos.size(); ++i){
            storePoint(a_pos[i], cl_pos[i], cd_pos[i]);
        }


        for (std::size_t i = 0; i < alphas_.size(); ++i)
            polar_ex_.pts.push_back( { alphas_[i], cls_[i], cds_[i], 0.0 } );
    }

    inline void storePoint(double a, double cl, double cd)
    {
        alphas_.push_back(a);
        cls_.push_back(cl);
        cds_.push_back(cd);
    }
};