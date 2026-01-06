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
        bool output_latest;
        bool output_np;
        double residual_limit{};
        bool run_state{};
        bool parallel_mode{};     // false for phy-parallel, true for dv-parallel
        int converge_state{};
        double mesh_scale{};
    public:
        explicit BasicSolver(ArgParser &parser);

        BasicSolver(ArgParser &parser, Config &config);

        virtual void update_config();

        bool get_run_state() const { return run_state; };

        int max_step();

        int write_interval();
    };

    /// Residual
    const int residual_interval = 100;

    Scalar residual(Field<Scalar> &_old, Field<Scalar> &_new);

    Vector residual(Field<Vector> &_old, Field<Vector> &_new);

    typedef std::vector<Field<Scalar>> DistributionFunction;

    /// limiter
    Scalar venkata_limiter(
            MESO::Field<Scalar> &f_field, MESO::Scalar dw, fvmMesh::Cell &cell,
            MESO::Scalar venkata_k, Scalar eps=1.0e-10
    );

/// add Solver here
#include "solver/boltzmann/dugks.h"
#include "solver/boltzmann/cdugks.h"
#include "solver/boltzmann/cdugks_shakhov.h"

/// handle Solver
    template<class SolverClass>
    int handle_solver(MESO::ArgParser &parser, MESO::Solver::Config &config);

    extern int solver_state;
}


#endif //MESO_SOLVER_H
