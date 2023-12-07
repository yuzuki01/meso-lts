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

Physical::MacroVars Physical::strvec_to_macro_vars(const string_vector &str_vec) {
    Physical::MacroVars result;
    if (str_vec.size() < 9) {
        std::string err_str = "strvec_to_macro_vars length error";
        warn_println(err_str);
        throw std::invalid_argument(err_str);
    }
    result.density = stod(str_vec[1]);
    result.temperature = stod(str_vec[2]);
    result.velocity = {stod(str_vec[3]), stod(str_vec[4]), stod(str_vec[5])};
    result.heat_flux = {stod(str_vec[6]), stod(str_vec[7]), stod(str_vec[8])};
    return result;
}
