#include "solver.h"

TP_func
LeastSquare generate_least_square(const int &cell_key, MESH::ListMesh &mesh) {
    LeastSquare result;
    const int dimension = mesh.dimension();
    auto &cell = mesh.get_cell(cell_key);
    const int near_cell_num = cell.near_cell_key.size();
    double Sxx, Sxy, Sxz, Syy, Syz, Szz;
    result.dr.resize(near_cell_num);
    result.weight.resize(near_cell_num);
    result.dr.shrink_to_fit();
    result.weight.shrink_to_fit();
    Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
    for (int i = 0; i < near_cell_num; i++) {
        auto &near_cell = mesh.get_cell(cell.near_cell_key[i]);
        Vec3D dr = near_cell.position - cell.position;
        double wi = 1.0 / (dr * dr);
        result.dr[i] = dr;
        result.weight[i] = wi;
        // sum
        Sxx += wi * dr.x * dr.x;
        Sxy += wi * dr.x * dr.y;
        Sxz += wi * dr.x * dr.z;
        Syy += wi * dr.y * dr.y;
        Syz += wi * dr.y * dr.z;
        Szz += wi * dr.z * dr.z;
    }
    if (dimension == 2) {
        Sxz = Syz = 0.0;
        Szz = 1.0;
    }
    Mat3D a{
        Sxx, Sxy, Sxz,
        Sxy, Syy, Syz,
        Sxz, Syz, Szz,
        };
    result.C = a.I();
    if (dimension == 2) {
        result.C._matrix[0][0] = 0.0;
        result.C._matrix[0][1] = 0.0;
        result.C._matrix[0][2] = 0.0;
    }
    return result;
}

TP_func
LeastSquare generate_least_square(const std::string &cell_key, MESH::MapMesh &mesh) {
    LeastSquare result;
    const int dimension = mesh.dimension();
    auto &cell = mesh.get_cell(cell_key);
    const int near_cell_num = cell.near_cell_key.size();
    double Sxx, Sxy, Sxz, Syy, Syz, Szz;
    result.dr.resize(near_cell_num);
    result.weight.resize(near_cell_num);
    result.dr.shrink_to_fit();
    result.weight.shrink_to_fit();
    Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
    for (int i = 0; i < near_cell_num; i++) {
        auto &near_cell = mesh.get_cell(cell.near_cell_key[i]);
        Vec3D dr = near_cell.position - cell.position;
        double wi = 1.0 / (dr * dr);
        result.dr[i] = dr;
        result.weight[i] = wi;
        // sum
        Sxx += wi * dr.x * dr.x;
        Sxy += wi * dr.x * dr.y;
        Sxz += wi * dr.x * dr.z;
        Syy += wi * dr.y * dr.y;
        Syz += wi * dr.y * dr.z;
        Szz += wi * dr.z * dr.z;
    }
    if (dimension == 2) {
        Sxz = Syz = 0.0;
        Szz = 1.0;
    }
    Mat3D a{
            Sxx, Sxy, Sxz,
            Sxy, Syy, Syz,
            Sxz, Syz, Szz,
    };
    result.C = a.I();
    if (dimension == 2) {
        result.C._matrix[0][0] = 0.0;
        result.C._matrix[0][1] = 0.0;
        result.C._matrix[0][2] = 0.0;
    }
    return result;
}


/// limiter ventaka
double ventaka_phi(const Ventaka &_limiter, double omega) {
    if (_limiter.delta_ij == 0.0) return 1.0;
    double a, aa, ab, bb;
    if (_limiter.delta_ij > 0.0) {
        a =_limiter.Wi_max - _limiter.Wi;
    } else {
        a = _limiter.Wi_min - _limiter.Wi;
    }
    aa = a * a;
    bb = _limiter.delta_ij * _limiter.delta_ij;
    ab = a * _limiter.delta_ij;
    return (aa + 2.0 * ab + omega) / (aa + 2.0 * bb + ab + omega);
}
