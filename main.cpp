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
    /// sample
    if (parser.parse_switch("test")) return test();

    if (debug_mode) debug_println("exit");

    return 0;
}

int test() {
    Mat3D a{3.0, 2.0, 8.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0};
    Mat3D b{1.0, 0.0, 0.0, 0.0, 2.0, 6.0, 1.0, 0.0, 3.0};
    std::cout << "a.det = " << a.det() << std::endl;
    std::cout << "Det(a.I * b) = " << (a.I() * b).det() << std::endl;
    auto it = a.I();
    it.info();
    return 0;
}
