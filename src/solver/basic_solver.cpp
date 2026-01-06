#include "solver/solver.h"

using namespace MESO::Solver;

BasicSolver::BasicSolver(ArgParser &parser) :
        parser(parser), config(parser.parse_param<std::string>("case", "<case-file>", true)) {
    case_name = config.get<std::string>("case-name", "unnamed", false);
    output_latest = config.get<bool>("output-latest", true, false);
    output_np = config.get<bool>("output-np", false, false);
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
    }
    residual_limit = config.get<double>("residual-limit", 1e-6, false);
    mesh_scale = config.get<double>("mesh-scale", 1.0, false);
}


BasicSolver::BasicSolver(MESO::ArgParser &parser, MESO::Solver::Config &config):
    parser(parser), config(config) {
    case_name = config.get<std::string>("case-name", "unnamed", false);
    output_latest = config.get<bool>("output-latest", true, false);
    output_np = config.get<bool>("output-np", false, false);
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
    }
    residual_limit = config.get<double>("residual-limit", 1e-6, false);
    mesh_scale = config.get<double>("mesh-scale", 1.0, false);
}

int BasicSolver::max_step() {
    return config.get<int>("max-step", 100000, false);
}

int BasicSolver::write_interval() {
    return config.get<int>("write-interval", 1000, false);
}

void BasicSolver::update_config() {
    config.update_config(false);
}
