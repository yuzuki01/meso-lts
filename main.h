#include "core.h"
#include "mesh.h"
#include "solver.h"

#include "sample.h"


template <class SOLVER>
int handle_solver(ConfigReader &config, ArgParser &parser)  {
    config.info();
    SOLVER solver(config, parser);
    solver.info();
    solver.init();
    solver.do_save();
    while (solver.continue_to_run) {
        solver.do_step();
    }
    solver.do_save();
    return 0;
}

/// explicit init
#ifdef SOLVER_DUGKS_INCOMPRESSIBLE
template int handle_solver<DUGKS_INCOMPRESSIBLE>(ConfigReader &config, ArgParser &parser);
#endif
