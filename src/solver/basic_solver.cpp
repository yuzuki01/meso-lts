#include "solver/solver.h"

using namespace MESO;

Solver::BasicSolver::BasicSolver(ArgParser &parser) :
        parser(parser), config(parser.parse_param<std::string>("case", "<case-file>", true)) {
    case_name = config.get<std::string>("case-name", "unnamed", false);
    output_np = config.get<bool>("output-np", false, false);
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
    }
}
