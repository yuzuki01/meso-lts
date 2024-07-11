#ifndef MESO_BOUNDARY_H
#define MESO_BOUNDARY_H

namespace MESO::Solver::Boundary {
    /// 不可压缩流动
    namespace Incompressible {
        Scalar solve_pressure(Scalar P0, Scalar Rho, Vector &U, Scalar limit=1e-13);
    }
    /// 可压缩流动
    namespace Compressible {
        Scalar solve_pressure(Scalar P0, Scalar R, Scalar T, Vector &U, Scalar gamma);
    }
}

#endif //MESO_BOUNDARY_H
