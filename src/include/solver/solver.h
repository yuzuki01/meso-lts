#ifndef MESO_SOLVER_H
#define MESO_SOLVER_H


#include "core/core.h"
#include "mesh/mesh.h"
#include "field/field.h"


namespace MESO::Solver {

    class BasicSolver {
    protected:
        FileIO::ParamReader config_;
        const String name_;
        Time time_;
        fvMesh mesh_;

    public:
        explicit BasicSolver(const String &filePath);

        /// Interfaces
        [[nodiscard]] const String &name() const;

        [[nodiscard]] const Time &time() const;

        [[nodiscard]] const fvMesh &mesh() const;

        /// Solution
        bool solution();

        /// File IO
        void output();
    };

    /// Solvers
#include "solver/boltzmann/cdugks.h"

}


#endif //MESO_SOLVER_H
