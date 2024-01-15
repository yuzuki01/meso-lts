#include "main.h"

int main(int argc, char** argv) {
    meso_init();

    ArgParser parser(argc, argv);
    debug_mode = parser.parse_switch("debug");
    if (parser.parse_switch("h") || parser.parse_switch("help")) return handle_help();

    std::string parsed_string;

    /// test
    if (parser.parse_switch("test")) return handle_test();

    /// mesh
    parsed_string = parser.parse_param<std::string>("parse_mesh", STRING_NULL);
    if (parsed_string != STRING_NULL) return handle_parse_mesh(parsed_string);

    /// case
    parsed_string = parser.parse_param<std::string>("case", STRING_NULL);
    if (parsed_string != STRING_NULL) {
        ConfigReader config(parsed_string);
#ifdef SOLVER_DUGKS_INCOMPRESSIBLE
        if (config.solver == "dugks@incompressible") return handle_solver<DUGKS_INCOMPRESSIBLE>(config, parser);
#endif
#ifdef SOLVER_DUGKS_SHAKHOV
        if (config.solver == "dugks@shakhov") return handle_solver<DUGKS_SHAKHOV>(config, parser);
#endif
    }

    if (debug_mode) debug_println("exit");

    return 0;
}
