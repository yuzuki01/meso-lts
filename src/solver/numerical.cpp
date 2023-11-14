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
    switch (dimension) {
        case 2: {
            double FM = Sxx * Syy - Sxy * Sxy;
            result.Cx = {Syy / FM, -Sxy / FM, 0.0};
            result.Cy = {-Sxy / FM, Sxx / FM, 0.0};
            result.Cz = {0.0, 0.0, 0.0};
            return result;
        }
        case 3: {
            double FM = -Sxz * Sxz * Syy + 2.0 * Sxy * Sxz * Syz - Sxx * Syz * Syz - Sxy * Sxy * Szz +
                        Sxx * Syy * Szz;
            double SxyyzNxzyy = Sxy * Syz - Sxz * Syy;
            double SxzyzNxyzz = Sxz * Syz - Sxy * Szz;
            double SyyzzNyzyz = Syy * Szz - Syz * Syz;
            double SxxzzNxzxz = Sxx * Szz - Sxz * Sxz;
            double SxyxzNxxyz = Sxy * Sxz - Sxx * Syz;
            double SxxyyNxyxy = Sxx * Syy - Sxy * Sxy;
            result.Cx = {
                    SyyzzNyzyz / FM, SxzyzNxyzz / FM, SxyyzNxzyy / FM
            };
            result.Cy = {
                    SxzyzNxyzz / FM, SxxzzNxzxz / FM, SxyxzNxxyz / FM
            };
            result.Cz = {
                    SxyyzNxzyy / FM, SxyxzNxxyz / FM, SxxyyNxyxy / FM
            };
            return result;
        }
        default:
            std::string s = "generate_least_square() caught invalid dimension.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
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
    switch (dimension) {
        case 2: {
            double FM = Sxx * Syy - Sxy * Sxy;
            result.Cx = {Syy / FM, -Sxy / FM, 0.0};
            result.Cy = {-Sxy / FM, Sxx / FM, 0.0};
            result.Cz = {0.0, 0.0, 0.0};
            return result;
        }
        case 3: {
            double FM = -Sxz * Sxz * Syy + 2.0 * Sxy * Sxz * Syz - Sxx * Syz * Syz - Sxy * Sxy * Szz +
                        Sxx * Syy * Szz;
            double SxyyzNxzyy = Sxy * Syz - Sxz * Syy;
            double SxzyzNxyzz = Sxz * Syz - Sxy * Szz;
            double SyyzzNyzyz = Syy * Szz - Syz * Syz;
            double SxxzzNxzxz = Sxx * Szz - Sxz * Sxz;
            double SxyxzNxxyz = Sxy * Sxz - Sxx * Syz;
            double SxxyyNxyxy = Sxx * Syy - Sxy * Sxy;
            result.Cx = {
                    SyyzzNyzyz / FM, SxzyzNxyzz / FM, SxyyzNxzyy / FM
            };
            result.Cy = {
                    SxzyzNxyzz / FM, SxxzzNxzxz / FM, SxyxzNxxyz / FM
            };
            result.Cz = {
                    SxyyzNxzyy / FM, SxyxzNxxyz / FM, SxxyyNxyxy / FM
            };
            return result;
        }
        default:
            std::string s = "generate_least_square() caught invalid dimension.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

