//
// Created by Lukas Campbell on 11/05/2025.
//

#pragma once

namespace xf {
// ─────────────────────────────────────────────────────────────────────────
//                    small utility funcs.
// ─────────────────────────────────────────────────────────────────────────
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


    // ─────────────────────────────────────────────────────────────────────────
    //                   config + runner for a single polar
    // ─────────────────────────────────────────────────────────────────────────
    struct Config {
        // --- xfoil config -----------------------------------------------------------
        bool headless = true;                      // disable plots in the container
        // --- airfoil ----------------------------------------------------------------
        std::string naca      = "2412";            // ignored if coord_file!=""
        std::string coord_file = "";               // rel. path
        // --- mesh -------------------------------------------------------------------
        int    panels      = 160;
        double te_le_ratio = 0.15;
        // --- flow -------------------------------------------------------------------
        double re      = 1e7;
        double mach    = 0.0;
        double n_crit  = 7.0;
        // --- solver -----------------------------------------------------------------
        int    ITER    = 300;
        double VACCEL  = 0.01;
        // --- α-sweep -----------------------------------------------------------------
        double alpha = 0.0;
        double alpha_neg = -6.0;
        double alpha_pos =  6.0;
        double d_alpha    =   1;
        // --- misc -------------------------------------------------------------------
        std::string tag    = "run";
        bool   save_cp      = false;
        static const std::unordered_map<std::string, std::string> save_paths;

        static std::string g_save_path(const std::string& key, const std::string& name = "") {
            if (save_paths.contains(key)) return get_project_root() + save_paths.at(key) + name;
            std::cerr << "[ERROR] Incorrect save path key " << key << std::endl;
            return "";
        }

        static std::string to_container_path(const std::string& host_path) {
            const std::string proj_root = get_project_root();
            const std::filesystem::path full(host_path);
            const std::filesystem::path rel = std::filesystem::relative(full, proj_root);
            return "/work/" + rel.generic_string();
        }

        // --- reproducibility ---------------------------------------------------------
        const std::string toJSON() const {
            std::ostringstream js;
            js << "{"
               << R"("naca":")" << naca << R"(",)"
               << R"("coord_file":")" << coord_file << R"(",)"
               << R"("panels":)" << panels << ","
               << R"("te_le_ratio":)" << te_le_ratio << ","
               << R"("Ncrit":)" << n_crit << ","
               << R"("ITER":)" << ITER << ","
               << R"("VACCEL":)" << VACCEL << ","
               << R"("alpha_neg":)" << alpha_neg << ","
               << R"("alpha_pos":)" << alpha_pos << ","
               << R"("d_alpha":)" << d_alpha
               << "}";
            return js.str();
        }
    };

    struct Point { double alpha, cl, cd, cm; };
    struct Polar { double Re = 0.0; std::vector<Point> pts; };
    // ─────────────────────────────────────────────────────────────────────────
    //                       docker-based runner
    // ─────────────────────────────────────────────────────────────────────────
    class XFOILRunner {
    public:
        /*
         *  image          : docker/podman image name (default pulls from hub)
         *  runtime        : docker or podman if ur on linux
         *  mountWholeCwd  : true → mount $PWD ; false → mount cfg.tmp_dir only
         */
        explicit XFOILRunner(std::string image = "thomaseizinger/xfoil",
                             std::string runtime = "docker",
                             bool mountWholeCwd = false)
            : _image(std::move(image)),
              _runtime(std::move(runtime)),
              _mountWholeCwd(mountWholeCwd)
        {
            for (const auto &val: Config::save_paths | std::views::values) {
                if (std::filesystem::path path = get_project_root() + val; !std::filesystem::exists(path)) {
                    std::cerr << "Creating directory: " << path << "\n";
                    std::filesystem::create_directories(path);
                }
            }
        }

    // ──────────────── routine helpers ───────────────────────────────────────────────────────
    // each function must always return to the XFOIL state (ie, terminal reads  XFOIL c>)

        static void toggle_headless(std::ofstream& ss) {
            ss << "PLOP" << "\n"
               << "  G"  << "\n"
               << "\n";
        }

        static void init_airfoil_geometry(std::ofstream& ss, const Config& cfg) {
            if(!cfg.coord_file.empty()) {
                ss << "LOAD " << cfg.coord_file << "\n";
            } else {
                ss << "NACA " << cfg.naca << "\n";
            }
        }

        static void configure_airfoil(std::ofstream& ss, const Config& cfg) {
            ss << "PPAR" << "\n"
               << "  N" << cfg.panels << "\n"
               << "  T" << cfg.te_le_ratio << "\n"
               << "\n"
               << "\n";

            ss << "PANE" << "\n";
        }

        static void configure_solver(std::ofstream& ss, const Config& cfg) {
            ss << "OPER" << "\n"
               << "  VISC" << std::scientific << cfg.re << "\n";

            if (cfg.mach > 0.0) {
                ss << "  MACH" << std::fixed << std::setprecision(4) << cfg.mach << "\n";
            }

            ss << "ITER" << cfg.ITER << "\n";

            ss << "VPAR" << "\n"
               << "  N " << cfg.n_crit << "\n"
               << "  VACCEL " << cfg.VACCEL << "\n"
               << "\n"
               << "\n";
        }

        static void set_alpha(std::ofstream& ss, const double& alpha) {
            ss << "OPER" << "\n"
               << "  ALFA" << std::fixed << std::setprecision(6) << alpha << "\n";
        }

        static void set_alpha_sweep(std::ofstream& ss, const Config& cfg) {
            for(double a = cfg.alpha_neg; a <= cfg.alpha_pos; a += cfg.d_alpha) {
                set_alpha(ss, a);
            }
        }

        static void quit(std::ofstream& ss) {
            ss << "\nQUIT\n";
            ss.close();
        }


        static void custom_command(std::ofstream& ss, const std::string& cmd_template, auto... args) {
            std::ostringstream formatted;
            formatted << std::vformat(cmd_template, std::make_format_args(args...));
            ss << formatted.str() << "\n";
        }


        static void save_airfoil(std::ofstream& ss, const Config& cfg, const std::string& save_cmd = "SAVE", const std::string& filename = "") {
            if (save_cmd.empty()) {
                std::cerr << "Error: No save command provided.\n";
                return;
            }

            if (save_cmd != "SAVE" && save_cmd != "PSAV" && save_cmd != "ISAV" && save_cmd != "MSAV") {
                std::cerr << "Error: Invalid save command. Use SAVE, PSAV, ISAV, or MSAV.\n";
                return;
            }

            const std::string name      = (filename.empty() ? cfg.naca : filename) + ".dat";
            const std::string hostPath  = Config::g_save_path("AIRFOIL", name);
            const std::string container = Config::to_container_path(hostPath);

            std::filesystem::create_directories(std::filesystem::path(hostPath).parent_path());

            std::cout << "airfoil filepath: " << container << "\n";
            ss << save_cmd << " \"" << container << "\"\n";
        }



    // ─────────── routine runners ─────────────────────────────────────────────────────────────
        bool run(const Config& cfg, Polar& polar, const std::function<void(std::ofstream&)>& routine, std::ostream* verb = nullptr) {
            std::filesystem::create_directories(Config::g_save_path("TMP"));
            const std::string base = Config::g_save_path("TMP") + cfg.tag;
            const std::string scr  = base + ".xfs";


            std::ofstream ss(scr);
            if(!ss) {
                std::cerr << "Cannot write " << scr << "\n";
                return false;
            }


            routine(ss);

            const std::string mount_path = std::filesystem::absolute(get_project_root()).generic_string();

            std::string cmd =
                _runtime + " run --rm -i " +
                "-v " + sh_quote(mount_path + ":/work") + " " +
                "-w /work " +
                _image  + " xfoil " +
                "< " + sh_quote(scr) + " 2>&1";

            if(verb) *verb << "[XFOIL] Running: " << cmd << "\n";

            std::string xfoil_output = shell(cmd, verb);

            if(verb) {
                *verb << "\n=== XFOIL Output ===\n";
                *verb << xfoil_output << "\n";
                *verb << "==================\n";
            }

            polar = parse_xfoil_output_to_polar(xfoil_output, cfg.re);
            // save_polar_file(polar, cfg, xfoil_output);

            std::filesystem::remove(scr);

            return !polar.pts.empty();
        }

    private:
        std::string _image, _runtime;
        bool _mountWholeCwd;

        static Polar parse_xfoil_output_to_polar(const std::string& output, double Re) {
            std::regex  rx_a_cl(R"(\s*a\s*=\s*([-+]?\d+\.?\d*)\s+CL\s*=\s*([-+]?\d+\.?\d*))",
                                std::regex::optimize);
            std::regex  rx_cm_cd(R"(\s*Cm\s*=\s*([-+]?\d+\.?\d*)\s+CD\s*=\s*([-+]?\d+\.?\d*))",
                                 std::regex::optimize);

            double cur_alpha = 0.0, cur_cl = 0.0;
            bool have_alpha = false;

            std::map<double, Point> latest;

            std::istringstream in(output);
            std::string line;

            while (std::getline(in, line)) {
                std::smatch m;
                if (std::regex_search(line, m, rx_a_cl)) {
                    cur_alpha   = std::stod(m[1]);
                    cur_cl      = std::stod(m[2]);
                    have_alpha  = true;
                    continue;
                }

                if (have_alpha && std::regex_search(line, m, rx_cm_cd)) {
                    Point p;
                    p.alpha = cur_alpha;
                    p.cl    = cur_cl;
                    p.cm    = std::stod(m[1]);
                    p.cd    = std::stod(m[2]);

                    latest[p.alpha] = p;
                    have_alpha = false;
                }
            }

            Polar out;
            out.Re = Re;
            out.pts.reserve(latest.size());
            for (auto& [_, pt] : latest) out.pts.push_back(pt);
            return out;
        }

        static void save_polar_file(const Polar& polar, const Config& cfg) {
            std::ostringstream filename;
            filename << Config::g_save_path("TMP");

            if (!cfg.coord_file.empty()) {
                std::filesystem::path p(cfg.coord_file);
                filename << p.stem().string();
            } else {
                filename << "naca" << cfg.naca;
            }

            filename << "_Re" << std::fixed << std::setprecision(0) << cfg.re;

            if (cfg.mach > 0.0) {
                filename << "_Mach" << std::fixed << std::setprecision(2) << cfg.mach;
            }

            filename << ".pol";

            std::ofstream out(filename.str());
            if (!out) {
                std::cerr << "Failed to create polar file: " << filename.str() << "\n";
                return;
            }

            out << "  \n";
            out << "       XFOIL         Version 6.99\n";
            out << "  \n";
            out << " Calculated polar for: ";
            if (!cfg.coord_file.empty()) {
                out << std::filesystem::path(cfg.coord_file).stem().string();
            } else {
                out << "NACA " << cfg.naca;
            }
            out << "                                       \n";
            out << "  \n";
            out << " 1 1 Reynolds number fixed          Mach number fixed         \n";
            out << "  \n";
            out << " xtrf =   1.000 (top)        1.000 (bottom)  \n";
            out << " Mach =   " << std::fixed << std::setprecision(3) << cfg.mach
                << "     Re =     " << std::scientific << std::setprecision(3) << cfg.re / 1e6
                << " e 6     Ncrit =   " << std::fixed << std::setprecision(3) << cfg.n_crit
                << "\n";
            out << "  \n";
            out << "   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr\n";
            out << "  ------ -------- --------- --------- -------- -------- --------\n";

            for (const auto& pt : polar.pts) {
                out << std::fixed << std::setw(8) << std::setprecision(3) << pt.alpha
                    << std::fixed << std::setw(9) << std::setprecision(4) << pt.cl
                    << std::fixed << std::setw(10) << std::setprecision(5) << pt.cd
                    << std::fixed << std::setw(10) << std::setprecision(5) << 0.0  // cp not available from output
                    << std::fixed << std::setw(9) << std::setprecision(4) << pt.cm
                    << std::fixed << std::setw(9) << std::setprecision(4) << 1.0   // top_xtr
                    << std::fixed << std::setw(9) << std::setprecision(4) << 1.0   // bot_xtr
                    << "\n";
            }

            out.close();
            std::cout << "Saved polar file: " << filename.str() << "\n";
        }
    };

    struct PolarFamily {
        std::vector<Polar> pol;

        bool save_csv(const Config& meta, const std::string& file_name = "") const {
            if(pol.empty()) {
                std::cerr << "Error: No polar data to save to CSV\n";
                return false;
            }

            std::filesystem::path filepath = meta.g_save_path("POLAR_TABLE") + (file_name.empty() ? (meta.naca + "_table.csv") : file_name);
            try {
                if (filepath.has_parent_path()) {
                    std::filesystem::create_directories(filepath.parent_path());
                }
            } catch (const std::exception& e) {
                std::cerr << "Failed to create directory for CSV file: " << e.what() << "\n";
                return false;
            }

            std::ofstream out(filepath);
            if(!out) {
                std::cerr << "Error: Failed to create file: " << filepath << "\n";
                return false;
            }

            out << "# " << meta.toJSON() << "\n";
            out << "# Re, alpha, Cl, Cd, Cm\n";
            for(const auto& p:pol) {
                for(const auto& q:p.pts) {
                    out << p.Re << "," << q.alpha << "," << q.cl << "," << q.cd << "," << q.cm << "\n";
                }
            }
            std::cout << "Successfully wrote " << file_name << " with data for " << pol.size() << " Reynolds numbers\n";
            return true;
        }
    };

    struct Viterna {
        static void extend(Polar& p) {
            if(p.pts.empty()) return;
            auto deg2rad = [](double a) {
                return a * M_PI / 180.0;
            };

            auto it = std::max_element(p.pts.begin(),p.pts.end(), [](auto&a, auto&b) {
                return fabs(a.cl) < fabs(b.cl);
            } );
            const double CLmax = it->cl;
            const double CD90  = 1.11 + 0.018*fabs(it->alpha);

            auto viterna = [&](double a_deg) {
                double a = deg2rad(a_deg);
                double s = sin(a), c = cos(a);
                double cl = CD90 * 0.5 * s * c + CLmax * (s * s) / c;
                double cd = CD90 * s * s + CLmax * fabs(s) * c;
                return std::pair<double,double>{cl,cd};
            };

            double max_alpha = p.pts.back().alpha;
            double min_alpha = p.pts.front().alpha;

            if (max_alpha < 180.0) {
                double step = 2.0;
                for(double a = max_alpha + step; a <= 180.0; a += step) {
                    auto [cl,cd] = viterna(a);
                    p.pts.push_back({a, cl, cd, 0.0});
                }
            }

            if (min_alpha > -180.0) {
                double step = 2.0;
                for(double a = min_alpha - step; a >= -180.0; a -= step) {
                    auto [cl,cd] = viterna(a);
                    p.pts.push_back({a, cl, cd, 0.0});
                }
            }

            std::sort(p.pts.begin(),p.pts.end(), [](auto&a, auto&b) {
                return a.alpha<b.alpha;
            });
        }
    };

    class AirfoilTable {
    public:
        void build(const PolarFamily& fam) {
            _Re.clear(); _α.clear(); _cl.clear(); _cd.clear();
            if(fam.pol.empty()) return;

            std::unordered_map<double,bool> alphamap;
            for(const auto& p:fam.pol) {
                for(const auto& pt:p.pts) {
                    alphamap[pt.alpha] = true;
                }
            }

            for(auto& kv:alphamap) {
                _α.push_back(kv.first);
            }

            std::sort(_α.begin(),_α.end());

            _Re.reserve(fam.pol.size());

            for(const auto& p:fam.pol) {
                _Re.push_back(p.Re);
            }

            std::sort(_Re.begin(),_Re.end());

            const size_t R=_Re.size(), A=_α.size();

            _cl.assign(R*A,0.0); _cd.assign(R*A,0.0);

            auto idx = [&](size_t i,size_t j) {
                return i*A+j;
            };

            for(size_t i = 0; i < R; ++i) {
                const auto& pol = fam.pol[i];
                for(const auto& pt:pol.pts) {
                    size_t j = std::lower_bound(_α.begin(), _α.end(), pt.alpha) - _α.begin();
                    _cl[idx(i,j)] = pt.cl;
                    _cd[idx(i,j)] = pt.cd;
                }
                for(size_t j=1;j<A;++j) {
                    if(_cl[idx(i,j)]==0.0) _cl[idx(i,j)] = _cl[idx(i,j-1)];
                    if(_cd[idx(i,j)]==0.0) _cd[idx(i,j)] = _cd[idx(i,j-1)];
                }
            }
        }
        double cl(double Re,double alpha) const { return interp(_cl,Re,alpha); }
        double cd(double Re,double alpha) const { return interp(_cd,Re,alpha); }

    private:
        std::vector<double> _Re,_α,_cl,_cd;
        double interp(const std::vector<double>& arr,double Re,double a) const {
            if(_Re.empty()) return 0.0;

            size_t i = std::lower_bound(_Re.begin(), _Re.end(), Re) - _Re.begin();
            if(i == _Re.size()) i = _Re.size() - 1;

            size_t j = std::lower_bound(_α.begin(), _α.end(), a) - _α.begin();
            if(j == _α.size()) j = _α.size() - 1;

            return arr[i * _α.size() + j];
        }
    };

    inline int demo() {
        XFOILRunner xf;
        Config cfg;

        cfg.naca = "2412";
        cfg.d_alpha = 1.0;
        cfg.alpha_neg = 0.0;
        cfg.alpha_pos = 10.0;
        std::vector<double> re_list = {5e4, 1e5, 2e5};

        PolarFamily fam;

        for(size_t k = 0; k < re_list.size(); ++k) {
            cfg.re  = re_list[k];
            cfg.tag = "naca2412_Re" + std::to_string(int(cfg.re));
            Polar p;

            /*if(!xf.run(cfg,&std::cout, p)) {
                std::cerr<<"XFOIL failed to generate points for Re=" << cfg.re << "\n";
                continue;
            }*/

            std::cout << "Generated " << p.pts.size() << " points for Re=" << cfg.re << "\n";

            // Viterna::extend(p);
            fam.pol.push_back(std::move(p));
        }

        AirfoilTable tbl;
        tbl.build(fam);


        if (!fam.save_csv(cfg)) {
            std::cerr << "Failed to save CSV file\n";
        }

        if (!fam.pol.empty() && !fam.pol[0].pts.empty()) {
            std::cout << "First polar has " << fam.pol[0].pts.size() << " points\n";
            std::cout << "Cl at Re="<<fam.pol[0].Re<<", α=0° = "<<tbl.cl(fam.pol[0].Re,0.0)<<"\n";
        }

        return 0;
    }

} // namespace xf