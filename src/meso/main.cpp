#include "meso.h"


int main(int argc, char **argv) {
    MESO::ArgParser parser(argc, argv);

    if (parser.parse_switch("debug")) {
        logger.level = -1;
    }
    logger.debug << "Running in debug mode." << std::endl;

    if (parser.parse_switch("h")) {
        MESO::MPI::Initialize(&argc, &argv);
        MESO::help();
        MESO::MPI::Finalize();
        return 0;
    }

    if (parser.parse_switch("v") || parser.parse_switch("version")) {
        MESO::version();
        return 0;
    }

    if (parser.parse_switch("mpi-test")) {
        MESO::MPI::Initialize(&argc, &argv);
        MESO::MPI::Finalize();
        return 0;
    }

    MESO::String parse_string;

    parse_string = parser.parse_param<MESO::String>("mesh", "<mesh-file>", false);
    if (parse_string != "<mesh-file>") {
        MESO::MPI::Initialize(&argc, &argv);
        auto code = MESO::handle_mesh(parse_string,
                                      parser.parse_param<double>("scale", 1.0, false));
        MESO::MPI::Finalize();
        return code;
    }

    parse_string = parser.parse_param<MESO::String>("case", "<case-file>", false);
    if (parse_string != "<case-file>") {
        auto solver = parser.parse_param<MESO::String>("solver", "<solver>", false);
        /// MPI init
        MESO::MPI::Initialize(&argc, &argv);
        int solver_status;
        MESO::Solver::Config config(parse_string);
        if (solver == "<solver>") solver = config.get<std::string>("solver", "<solver>", false);
        if (solver == "dugks@incompressible")
            solver_status = MESO::Solver::handle_solver<MESO::Solver::DUGKS>(parser, config);
        else if (solver == "cdugks@incompressible")
            solver_status = MESO::Solver::handle_solver<MESO::Solver::CDUGKS>(parser, config);
        else if (solver == "cdugks@shakhov")
            solver_status = MESO::Solver::handle_solver<MESO::Solver::CDUGKS_SHAKHOV>(parser, config);
        else {
            logger.warn << "Load solver: " << solver << " failed." << std::endl;
            solver_status = 0;
        }
        MESO::MPI::Finalize();
        return solver_status;
    }

    /// run with no params
    logger.warn << "type \"./meso-mpi -h\" to get help." << std::endl;
    return 0;
}
