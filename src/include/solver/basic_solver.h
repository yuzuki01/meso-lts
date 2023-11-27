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
