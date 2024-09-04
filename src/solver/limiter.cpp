#include "solver/solver.h"


using namespace MESO;


Scalar Solver::venkata_limiter(Field<MESO::Scalar> &f_field, MESO::Scalar df,
                               fvmMesh::Face &face, fvmMesh::Cell &cell,
                               MESO::Scalar venkata_k) {
    /// find f_min and f_max
    auto f = f_field[cell.id],
            f_min = f_field[cell.id],
            f_max = f_field[cell.id];
    for (auto neighbor_id: cell.neighbors) {
        auto f_i = f_field[neighbor_id];
        if (f_i < f_min) f_min = f_i;
        if (f_i > f_max) f_max = f_i;
    }
    auto &mesh = *f_field.get_mesh();
    auto kh = venkata_k * pow(cell.volume, 1.0 / mesh.dimension());
    auto omega = kh * kh * kh;
    if (df == 0.0) return 1.0;
    Scalar a, aa, ab, bb;
    if (df > 0.0) {
        a = f_max - f;
    } else {
        a = f_min - f;
    }
    aa = a * a;
    ab = a * df;
    bb = df * df;
    return (aa + 2.0 * ab + omega) / (aa + 2.0 * bb + ab + omega);
}
