#include "main.h"

int main(int argc, char** argv) {
    meso_init();

    ArgParser parser(argc, argv);
    debug_mode = parser.parse_switch("debug");
    if (parser.parse_switch("h") || parser.parse_switch("help")) return handle_help();

    std::string parsed_string;

    /// mesh
    parsed_string = parser.parse_param<std::string>("parse_mesh", STRING_NULL);
    if (parsed_string != STRING_NULL) return handle_parse_mesh(parsed_string);

    /// case
    parsed_string = parser.parse_param<std::string>("case", STRING_NULL);
    if (parsed_string != STRING_NULL) {
        ConfigReader config(parsed_string);
        if (config.solver == "dugks@incompressible") return handle_solver<DUGKS_INCOMPRESSIBLE>(config, parser);
        else if (config.solver == "cdugks@incompressible") return handle_solver<CDUGKS_INCOMPRESSIBLE>(config, parser);
        else if (config.solver == "dugks@shakhov") return handle_solver<DUGKS_SHAKHOV>(config, parser);
    }

    if (debug_mode) debug_println("exit");

    return 0;
}
