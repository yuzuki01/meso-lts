#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    Solver::BasicSolver solver("system/config");

    while (solver.solution()) {
        logger.info << solver.name() << " - Time: " << solver.time().name() << std::endl;
    }

#include "mesoFinalize.h"

    return 0;
}
