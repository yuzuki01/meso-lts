#ifndef MESO_SOLVER_H
#define MESO_SOLVER_H

#include "mesh/mesh.h"

#include "solver/config.h"


namespace MESO::Solver {

    class BasicSolver {
    protected:
        ArgParser &parser;
        Config config;
    public:
        explicit BasicSolver(ArgParser &parser) :
                parser(parser), config(parser.parse_param<std::string>("case", "<case-file>", true)) {};
    };

/// add Solver here

    typedef std::vector<Field<Scalar>> DistributionFunction;

#include "solver/boltzmann/cdugks.h"

}


#endif //MESO_SOLVER_H
