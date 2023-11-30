template <class SOLVER>
int handle_solver(ConfigReader &config, ArgParser &parser)  {
    config.info();
    SOLVER solver(config, parser);
    solver.init();
    if (!solver.continue_to_run) {
        warn_println("Solver init failed.");
        return -1;
    }
    solver.info();
    solver.do_save();
    while (solver.continue_to_run) {
        solver.do_step();
    }
    solver.do_save();
    return 0;
}

/// explicit init
#ifdef SOLVER_DUGKS_INCOMPRESSIBLE
template int handle_solver<DUGKS_INCOMPRESSIBLE>(ConfigReader &config, ArgParser &parser);
#endif
#ifdef SOLVER_DUGKS_SHAKHOV
template int handle_solver<DUGKS_SHAKHOV>(ConfigReader &config, ArgParser &parser);
#endif

/// parse mesh
int handle_parse_mesh(const std::string &path);

/// help
int handle_help();
