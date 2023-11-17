#include "main.h"

int main(int argc, char** argv) {
    meso_init();

    ArgParser parser(argc, argv);
    debug_mode = parser.parse_switch("debug");

    std::string parsed_string;
    /// case
    parsed_string = parser.parse_param<std::string>("case", STRING_NULL);
    if (parsed_string != STRING_NULL) {
        ConfigReader config("./config/" + parsed_string);
        if (config.solver == "dugks@incompressible") return handle_solver<DUGKS_INCOMPRESSIBLE>(config, parser);
    }

    if (debug_mode) debug_println("exit");

    return 0;
}
