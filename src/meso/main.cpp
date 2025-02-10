#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    Solver::BasicSolver solver(
            FileIO::ParamReader("system/config")
            );

    while (solver.solution()) {
        logger.info << solver.name() << " - Time: " << solver.time().name() << std::endl;
    }

    MPI::Finalize();

    /// run with no params
    logger.warn << "type \"./meso-mpi -help\" to get help." << std::endl;
    return 0;
}
