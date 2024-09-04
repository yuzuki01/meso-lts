#ifndef MESO_SOLVER_H
#define MESO_SOLVER_H

#include "mesh/mesh.h"

#include "solver/config.h"

#include "solver/boundary.h"

/// DVS
#include "solver/dvs/dvs_reader.h"

namespace MESO::MPI {
    MPI_TaskObject DVS_partition(fvmMesh::Mesh &mesh);
}


namespace MESO::Solver {

    fvmMesh::Mesh generate_newton_cotes(int dimension, int n, int mount, double scale);

    class BasicSolver {
    protected:
        ArgParser &parser;
        Config config;
        std::string case_name;
        bool output_np;
        double residual_limit{};
        bool run_state{};
        int converge_state{};
        double mesh_scale{};
    public:
        explicit BasicSolver(ArgParser &parser);

        BasicSolver(ArgParser &parser, Config &config);

        bool get_run_state() const { return run_state; };
    };

    /// Residual
    const int residual_interval = 100;

    Scalar residual(Field<Scalar> &_old, Field<Scalar> &_new);
    Vector residual(Field<Vector> &_old, Field<Vector> &_new);

    typedef std::vector<Field<Scalar>> DistributionFunction;

    /// limiter
    Scalar venkata_limiter(Field<Scalar> &f_field, Scalar df, fvmMesh::Face &face, fvmMesh::Cell &cell, Scalar venkata_k);

/// add Solver here
#include "solver/boltzmann/cdugks.h"
#include "solver/boltzmann/cdugks_shakhov.h"

/// handle Solver
    template <class SolverClass>
    int handle_solver(MESO::ArgParser &parser, MESO::Solver::Config &config);

    extern int solver_state;

    template <class SolverClass>
    void solver_interrupt(int signum);
}


#endif //MESO_SOLVER_H
