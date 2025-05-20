//
// Created by Lukas Campbell on 11/05/2025.
//

#pragma once

namespace xf {
    using namespace units::angle;
    using namespace units::dimensionless;

    struct Config {
        // --- xfoil config -----------------------------------------------------------
        bool headless = true;                      // disable plots in the container
        // --- airfoil ----------------------------------------------------------------
        std::string naca       = "2412";           // ignored if coord_file!=""
        std::string coord_file;                    // rel. path
        // --- mesh -------------------------------------------------------------------
        int    panels      = 200;
        double te_le_ratio = 0.15;
        // --- flow -------------------------------------------------------------------
        dimensionless_t re     = 1e7;
        dimensionless_t mach   = 0.0;
        double n_crit = 9.0;
        // --- solver -----------------------------------------------------------------
        int    ITER   = 500;
        double VACCEL = 0.01;
        // --- α-sweep -----------------------------------------------------------------
        degree_t alpha = degree_t(0.0);
        degree_t alpha_neg = degree_t(-6.0);
        degree_t alpha_pos = degree_t(8.0);
        degree_t d_alpha = degree_t(0.25);
        // --- misc -------------------------------------------------------------------
        bool        save_cp  = false;
        std::string tag      = "run";

        static const std::unordered_map<std::string, std::string> save_paths;
        // ----------------------------------------------------------------------------



        static std::string g_save_path(const std::string& key, const std::string& name = "") {
            if (save_paths.contains(key)) return Helpers::get_project_root() + save_paths.at(key) + name;
            std::cerr << "[ERROR] Incorrect save path key " << key << std::endl;
            return "";
        }

        static std::string to_container_path(const std::string& host_path) {
            const std::string proj_root = Helpers::get_project_root();
            const std::filesystem::path full(host_path);
            const std::filesystem::path rel = std::filesystem::relative(full, proj_root);
            return "/work/" + rel.generic_string();
        }

        // --- reproducibility ---------------------------------------------------------
        std::string toJSON() const {
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




    struct Point {
        degree_t alpha;
        dimensionless_t c_l;
        dimensionless_t c_d;
        dimensionless_t c_m;
    };

    struct Polar {
        dimensionless_t re = 0.0;
        std::vector<Point> pts;
    };


    class XFOILRunner {
    public:
        explicit XFOILRunner(std::string image = "thomaseizinger/xfoil", std::string runtime = "auto")
            : image(std::move(image))
        {
            if (runtime == "auto") {
                #ifdef _WIN32
                    this->runtime = "cmd";
                #elif __APPLE__
                                this->runtime = "docker";
                #elif __linux__
                                this->runtime = "docker";
                #else
                                std::cerr << "Warning: Unknown OS, defaulting to docker runtime\n";
                                this->runtime = "docker";
                #endif
            } else {
                this->runtime = runtime;
            }

            for (const auto &val: Config::save_paths | std::views::values) {
                if (std::filesystem::path path = Helpers::get_project_root() + val; !std::filesystem::exists(path)) {
                    std::cerr << "Creating directory: " << path << "\n";
                    std::filesystem::create_directories(path);
                }
            }
        }
        // ──────────────── install/setup helpers ───────────────────────────────────────────────────────

        std::string get_runtime() const {
            return runtime;
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
            ss << "\n";
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
               << "  VISC" << std::scientific << cfg.re.value() << "\n";

            if (cfg.mach > 0.0) {
                ss << "  MACH" << std::fixed << std::setprecision(4) << cfg.mach.value() << "\n";
            }

            ss << "ITER" << cfg.ITER << "\n";

            ss << "VPAR" << "\n"
               << "  N " << cfg.n_crit << "\n"
               << "  VACCEL " << cfg.VACCEL << "\n"
               << "\n"
               << "\n";
        }

        static void dump_coordinates(std::ofstream& ss) {
            ss << "PANE\n";
            ss << "DUMP\n";
        }


        static void set_alpha(std::ofstream& ss, const degree_t& alpha) {
            ss << "OPER" << "\n"
               << "  ALFA" << std::fixed << std::setprecision(6) << alpha.value() << "\n"
               << "\n";
        }

        static void set_alpha_sweep(std::ofstream& ss, const Config& cfg) {
            for(degree_t a = cfg.alpha_neg; a <= cfg.alpha_pos; a += cfg.d_alpha) {
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
            const std::string name = (filename.empty() ? cfg.naca : filename) + ".dat";
            const std::string hostPath = Config::g_save_path("AIRFOIL", name);
            const std::string container = Config::to_container_path(hostPath);

            ss << save_cmd << " \"" << container << "\"\n";
            ss << "\n";
        }





    // ─────────── routine runners ─────────────────────────────────────────────────────────────
        bool run(const Config& cfg, Polar& polar, const std::function<void(std::ofstream&)>& routine, std::ostream* verb = nullptr) const {
            std::filesystem::create_directories(Config::g_save_path("TMP"));
            const std::string base = Config::g_save_path("TMP") + cfg.tag;
            const std::string scr  = base + ".xfs";

            std::ofstream ss(scr);
            if(!ss) {
                std::cerr << "Cannot write " << scr << "\n";
                return false;
            }

            routine(ss);

            std::string cmd;
            if (runtime == "docker") {
                const std::string mount_path = std::filesystem::absolute(Helpers::get_project_root()).generic_string();
                std::string container_id = "xfoil_run_" + std::to_string(std::rand());

                cmd = runtime + " run --name " + container_id + " --rm -i " +
                      "-v " + Helpers::sh_quote(mount_path + ":/work:rw") + " " +
                      "-w /work " +
                      image + " xfoil " +
                      "< " + Helpers::sh_quote(scr) + " 2>&1";


                // const std::string mount_path = std::filesystem::absolute(Helpers::get_project_root()).generic_string();
                //
                // cmd = runtime + " run --rm -i " +
                //       "-v " + Helpers::sh_quote(mount_path + ":/work:rw") + " " +
                //       "-w /work " +
                //       image + " xfoil " +
                //       "< " + Helpers::sh_quote(scr) + " 2>&1";
            } else if (runtime == "cmd") {
                std::string script_content = adjust_paths_for_local_execution(scr);
                if (!script_content.empty()) {
                    std::ofstream new_script(scr, std::ios::trunc);
                    new_script << script_content;
                    new_script.close();
                }

                std::string xfoil_exe_path = Helpers::get_project_root() + "/xfoil/xfoil.exe";
                std::string xfoil_cmd = xfoil_exe_path.empty() ? "xfoil.exe" : xfoil_exe_path;
                cmd = xfoil_cmd + " < " + scr + " 2>&1";
            } else {
                std::cerr << "Unsupported runtime: " << runtime << "\n";
                return false;
            }

            if(verb) *verb << "[XFOIL] Running: " << cmd << "\n";

            std::string xfoil_output = Helpers::shell(cmd, verb);

            if(verb) {
                *verb << "\n=== XFOIL Output ===\n";
                *verb << xfoil_output << "\n";
                *verb << "==================\n";
            }

            polar = parse_xfoil_output_to_polar(xfoil_output, cfg.re);
            std::filesystem::remove(scr);

            return !polar.pts.empty();
        }
    private:
        std::string image, runtime;

        static Polar parse_xfoil_output_to_polar(const std::string& output, double re) {
            std::regex rx_a_cl(R"(\s*a\s*=\s*([-+]?\d+\.?\d*)\s+CL\s*=\s*([-+]?\d+\.?\d*))",
                               std::regex::optimize);
            std::regex rx_cm_cd(R"(\s*Cm\s*=\s*([-+]?\d+\.?\d*)\s+CD\s*=\s*([-+]?\d+\.?\d*))",
                                std::regex::optimize);

            double cur_alpha = 0.0, cur_cl = 0.0;
            bool have_alpha = false;

            std::map<degree_t, Point> latest;

            std::istringstream in(output);
            std::string line;

            while (std::getline(in, line)) {
                std::smatch m;
                if (std::regex_search(line, m, rx_a_cl)) {
                    cur_alpha = std::stod(m[1]);
                    cur_cl = std::stod(m[2]);
                    have_alpha = true;
                    continue;
                }

                if (have_alpha && std::regex_search(line, m, rx_cm_cd)) {
                    Point p;
                    p.alpha = degree_t(cur_alpha);
                    p.c_l = dimensionless_t(cur_cl);
                    p.c_m = dimensionless_t(std::stod(m[1]));
                    p.c_d = dimensionless_t(std::stod(m[2]));

                    latest[p.alpha] = p;
                    have_alpha = false;
                }
            }

            Polar out;
            out.re = dimensionless_t(re);
            out.pts.reserve(latest.size());
            for (auto &pt: latest | std::views::values) out.pts.push_back(pt);
            return out;
        }

        static std::string adjust_paths_for_local_execution(const std::string& script_path) {
            std::ifstream input(script_path);
            if (!input) {
                std::cerr << "Cannot read " << script_path << " for path adjustment\n";
                return "";
            }

            std::string content, line;
            while (std::getline(input, line)) {
                size_t pos = line.find("/work/");
                if (pos != std::string::npos) {
                    std::string project_root = Helpers::get_project_root();
                    line.replace(pos, 6, project_root);

                    #ifdef _WIN32
                        std::replace(line.begin(), line.end(), '/', '\\');
                    #endif
                }
                content += line + "\n";
            }
            input.close();
            return content;
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

            filename << "_Re" << std::fixed << std::setprecision(0) << cfg.re.value();

            if (cfg.mach > 0.0) {
                filename << "_Mach" << std::fixed << std::setprecision(2) << cfg.mach.value();
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
            out << " Mach =   " << std::fixed << std::setprecision(3) << cfg.mach.value()
                << "     Re =     " << std::scientific << std::setprecision(3) << cfg.re.value() / 1e6
                << " e 6     Ncrit =   " << std::fixed << std::setprecision(3) << cfg.n_crit
                << "\n";
            out << "  \n";
            out << "   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr\n";
            out << "  ------ -------- --------- --------- -------- -------- --------\n";

            for (const auto&[alpha, c_l, c_d, c_m] : polar.pts) {
                out << std::fixed << std::setw(8) << std::setprecision(3) << alpha
                    << std::fixed << std::setw(9) << std::setprecision(4) << c_l.value()
                    << std::fixed << std::setw(10) << std::setprecision(5) << c_d.value()
                    << std::fixed << std::setw(10) << std::setprecision(5) << 0.0  // cp not available from output
                    << std::fixed << std::setw(9) << std::setprecision(4) << c_m.value()
                    << std::fixed << std::setw(9) << std::setprecision(4) << 1.0   // top_xtr
                    << std::fixed << std::setw(9) << std::setprecision(4) << 1.0   // bot_xtr
                    << "\n";
            }

            out.close();
            std::cout << "Saved polar file: " << filename.str() << "\n";
        }
    };
} // namespace xf