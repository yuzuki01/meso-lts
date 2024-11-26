#include "solver/solver.h"


using namespace MESO;


Scalar Solver::venkata_limiter(MESO::Field<Scalar> &f_field, MESO::Scalar dw,
                               fvmMesh::Cell &cell, MESO::Scalar venkata_k, Scalar eps) {
    if (fabs(dw) <= eps) return 1.0;
    auto &mesh = *f_field.get_mesh();
    auto kh = venkata_k * pow(cell.volume, 1.0 / static_cast<Scalar>(mesh.dimension()));
    auto omega = kh * kh * kh;
    auto w = f_field[cell.id];
    auto w_max = w;
    auto w_min = w;
    for (auto &neighbor: cell.neighbors) {
        auto w_n = f_field[neighbor];
        if (w_n > w_max) w_max = w_n;
        if (w_n < w_min) w_min = w_n;
    }
    Scalar a, aa, ab, bb;
    if (dw > 0.0) {
        a = w_max - w;
    } else {
        a = w_min - w;
    }
    aa = a * a;
    ab = a * dw;
    bb = dw * dw;
    return (aa + 2.0 * ab + omega) / (aa + 2.0 * bb + ab + omega);
}
