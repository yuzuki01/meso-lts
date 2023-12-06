/**
 * included by solver.h
 */

/// Basic Solver Class
class BasicSolver {
public:
    Logger logger;
    ConfigReader &config;
    ArgParser &parser;
    /// Solver Params
    int step{}, max_step{}, save_interval{};
    bool is_crashed{}, continue_to_run{}, stop_at_specific_time{};
    double simulate_time{}, stop_time{};

    explicit BasicSolver(ConfigReader &_config, ArgParser &_parser);
    void init();
    void info();
    void do_save();
    void do_step();
};

/// Physical Var
namespace PhysicalVar {
    struct MacroVars {
        /// conserved
        double density = 0.0, old_density = 0.0;
        double energy = 0.0, old_energy = 0.0;
        Vec3D momentum{0.0, 0.0, 0.0};
        /// primitive
        double temperature = 0.0;
        Vec3D velocity{0.0, 0.0, 0.0};
        Vec3D heat_flux{0.0, 0.0, 0.0};
    };
}
