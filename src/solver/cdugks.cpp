#include "solver/solver.h"

using namespace MESO;
using namespace Solver;


CDUGKS::CDUGKS(const MESO::String &filePath)
        : BasicSolver(filePath),
          DV_(config_.get<String>("DVMProperties",
                                  "mesh",
                                  "constant/dvs.neu"),
              time_,
              true),
          rho(mesh_),
          U(mesh_),
          gVol_(DV_.cellNum(), volScalarField(mesh_)),
          gSurf_(DV_.cellNum(), surfScalarField(mesh_))
          {
    mesh_.info();
}

/// Initialize
void CDUGKS::initialize() {

}

/// Solution
bool CDUGKS::solution() {
    time_.runStep();
    return time_.isStoppable();
}
