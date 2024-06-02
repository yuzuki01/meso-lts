#include "meso.h"


int main(int argc, char **argv) {
    /// MPI RUN
    MESO::MPI::Initialize(&argc, &argv);

    MESO::ArgParser parser(argc, argv);
    if (parser.parse_switch("debug")) {
        logger.level = -1;
        logger.debug << "Running in debug mode." << std::endl;
    }

    std::string parse_string;
    parse_string = parser.parse_param<std::string>("mesh", "<mesh-file>", false);
    if (parse_string != "<mesh-file>") {
        auto mesh = MESO::Mesh::load_gambit(parse_string);
        mesh.build_geom();
        mesh.info();
        mesh.output(parse_string);
        return 0;
    }

    parse_string = parser.parse_param<std::string>("case", "<case-file>", false);
    if (parse_string != "<case-file>") {
        omp_set_num_threads(parser.parse_param<int>("parallel", omp_get_max_threads(), true));
        int save_interval = parser.parse_param<int>("save-interval", 1000, false);
        MESO::Solver::CDUGKS solver(parser);
        solver.initial();
        solver.output();
        if (parser.parse_switch("not-run")) {
            MESO::MPI::Finalize();
            return 0;
        }
        for (int i = 0; i < parser.parse_param("max-step", 10000, false); ++i) {
            solver.do_step();
            if (solver.is_crashed) break;
            if (solver.step % save_interval == 0) solver.output();
        }
        solver.output();
        return 0;
    }
    MESO::MPI::Finalize();
    return 0;
}
