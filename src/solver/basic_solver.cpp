#include "solver.h"


BasicSolver::BasicSolver(ConfigReader &_config, ArgParser &_parser) : config(_config), parser(_parser) {
    /// Logger
    logger(config.solver);
    /// File
    create_dir("./result/" + config.name);
    /// Parser
    max_step = parser.parse_param<int>("max_step", 1e5);
    save_interval = parser.parse_param<int>("save_interval", 1000);
    step = parser.parse_param<int>("continue_step", 0);
    is_crashed = false;
}
