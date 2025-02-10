#include "solver/solver.h"


namespaceMESO


Solver::BasicSolver::BasicSolver(FileIO::ParamReader config)
        : config_(std::move(config)),
          name_(config_.get<String>("solver", "name", "EmptyName")),
          time_(Time(config_)),
          mesh_(fvMesh(
                  FileIO::BasicReader(config_.get<String>("solver",
                                                          "mesh",
                                                          "constant/mesh.neu")),
                  time_, true)) {

}


const String &Solver::BasicSolver::name() const {
    return name_;
}

const Time &Solver::BasicSolver::time() const {
    return time_;
}

const fvMesh &Solver::BasicSolver::mesh() const {
    return mesh_;
}

void Solver::BasicSolver::output() {
    logger.warn << "MESO::BasicSolver::output() was called: output grid file only" << std::endl;
    mesh_.output();
}

bool Solver::BasicSolver::solution() {
    logger.warn << "MESO::BasicSolver::solution() was called: nothing to be solved" << std::endl;
    time_.runStep();
    return time_.isStoppable();
}
