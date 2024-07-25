#ifndef MESO_BOUNDARY_H
#define MESO_BOUNDARY_H

namespace MESO::Solver::Boundary {
    /// 滑移壁面
    namespace SlipWall {
        Vector solve_slip_velocity(Vector &u_w, Vector &u_c,
                                   Vector &position_w, Vector &position_c,
                                   Vector &nv_interior, Scalar alpha, Scalar mfp);

        Scalar solve_slip_temperature(Scalar T_w, Scalar T_c,
                                      Vector &position_w, Vector &position_c,
                                      Vector &nv_interior, Scalar alpha, Scalar mfp,
                                      Scalar gamma, Scalar Pr);
    }
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
