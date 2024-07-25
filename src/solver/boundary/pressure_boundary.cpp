#include "solver/solver.h"


using namespace MESO::Solver::Boundary;

MESO::Scalar Incompressible::solve_pressure(Scalar P0, Scalar Rho, Vector &U, Scalar limit) {
    auto P = P0 - 0.5 * Rho * (U * U);
    if (P < limit) {
        std::cout << "Incompressible pressure boundary error: reached the limit: P=" << limit << std::endl;
        return limit;
    }
    return P;
}

MESO::Scalar Compressible::solve_pressure(Scalar P0, Scalar R, Scalar T, Vector &U, Scalar gamma) {
    auto Ma = U.magnitude() / sqrt(gamma * R * T);
    return P0 / pow(1.0 + (gamma - 1.0) / 2.0 * (Ma * Ma), gamma / (gamma - 1.0));
}
