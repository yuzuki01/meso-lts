#include "meso.h"


int main(int argc, char **argv) {
    MESO::ArgParser parser(argc, argv);
    
    if (parser.parse_switch("debug")) {
        logger.level = -1;
        logger.debug << "Running in debug mode." << std::endl;
    }

    MESO::String parse_string;
    parse_string = parser.parse_param<MESO::String>("mesh", "<mesh-file>", false);
    if (parse_string != "<mesh-file>") {
        auto mesh = MESO::Mesh::load_gambit(parse_string);
        mesh.build_geom();
        mesh.info();
        mesh.output(parse_string);
        return 0;
    }

    parse_string = parser.parse_param<MESO::String>("case", "<case-file>", false);
    if (parse_string != "<case-file>") {
        auto solver = parser.parse_param<MESO::String>("solver", "cdugks", false);
        if (solver == "cdugks") return MESO::Solver::handle_solver<MESO::Solver::CDUGKS>(&argc, &argv, parser);
        else return 0;
    }
    return 0;
}
