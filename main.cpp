#include "main.h"

int main(int argc, char** argv) {
    meso_init();

    ArgParser parser(argc, argv);
    debug_mode = parser.parse_switch("debug");

    if (debug_mode) {
        ConfigReader config("./config/demo.txt");
        config.info();
        DUGKS_INCOMPRESSIBLE solver(config, parser);
        // solver.info();
        solver.init();
        solver.do_save();
        while (solver.continue_to_run) {
            solver.do_step();
        }
        solver.do_save();
    }

    if (debug_mode) debug_println("exit");
    return 0;
}
