/**************************************
 * meso  - SOLVER module              *
 *  contains:                         *
 *         numerical, dugks solvers   *
 *               Nov 11, 2023  by MYC *
 **************************************/

#ifndef HEADER_SOLVER
#define HEADER_SOLVER

#ifndef HEADER_CORE
#include "core.h"
#endif

#ifndef HEADER_MESH
#include "mesh.h"
#endif

/// config reader
#include "solver/config_reader.h"
/// numerical
#include "solver/numerical.h"
/// Solver
#include "solver/basic_solver.h"    // basic

#include "solver/dugks/dugks_incompressible.h"


#endif  // HEADER_SOLVER
