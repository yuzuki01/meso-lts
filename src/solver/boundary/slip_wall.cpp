#include "solver/solver.h"


using namespace MESO::Solver::Boundary;


MESO::Vector SlipWall::solve_slip_velocity(MESO::Vector &u_w, MESO::Vector &u_c, MESO::Vector &position_w,
                                           MESO::Vector &position_c, MESO::Vector &nv_interior,
                                           MESO::Scalar alpha, MESO::Scalar mfp) {
    auto du_dn = (u_w - u_c) / ((position_w - position_c) * nv_interior);   // du_dn
    du_dn = du_dn - (du_dn * nv_interior) * nv_interior;    // (du_dn)_wall
    return u_w - 0.998 * ((2.0 - alpha) / alpha) * mfp * du_dn;
}

MESO::Scalar SlipWall::solve_slip_temperature(MESO::Scalar T_w, MESO::Scalar T_c, MESO::Vector &position_w,
                                              MESO::Vector &position_c, MESO::Vector &nv_interior, MESO::Scalar alpha,
                                              MESO::Scalar mfp, MESO::Scalar gamma, MESO::Scalar Pr) {
    auto dT_dn = (T_w - T_c) / ((position_w - position_c) * nv_interior);   // dT_dn
    return T_w - 0.982 * ((2 * gamma) / (gamma + 1)) * (mfp / Pr) * ((2.0 - alpha) / alpha) * dT_dn;
}
