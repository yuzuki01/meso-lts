#ifndef MESO_SOLVER_H
#define MESO_SOLVER_H

#include "mesh/mesh.h"

#include "solver/config.h"


namespace MESO::MPI {
    MPI_TaskObject DVS_partition(Mesh::Mesh &mesh);
}


namespace MESO::Solver {

    Mesh::Mesh generate_gauss_hermite(int dimension, int _n, double RT);

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
