#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    Solver::CDUGKS solver("system/config");

    solver.initialize();
    while (solver.solution()) {

    }

#include "mesoFinalize.h"

    return 0;
}
