#include "solver/solver.h"


typedef std::unordered_map<int, std::vector<double>> ValueListMap;

ValueListMap nc_weight = {
        {1, {1. / 2.,    1. / 2.}},
        {2, {1. / 6.,    4. / 6.,   1. / 6.}},
        {3, {1. / 8.,    3. / 8.,   3. / 8.,    1. / 8.}},
        {4, {7. / 90.,   16. / 45., 2. / 15.,   16. / 45.,  7. / 90.}},
        {5, {19. / 288., 25. / 96., 25. / 144., 25. / 144., 25. / 96., 19. / 288.}},
        {6, {41. / 840., 9. / 35.,  9. / 280.,  34. / 105., 9. / 280., 9. / 35., 41. / 840.}}
};

using namespace MESO;


Mesh::Mesh Solver::generate_newton_cotes(int dimension, int n, int mount, double scale) {
    int total_point = n * mount + 1;
    auto &p_weight = nc_weight[n];
    double dh = 2.0 * scale / mount;
    ScalarList coordinate(total_point), weight(total_point);
    Mesh::Mesh mesh;
    for (int i = 0; i < total_point; i++) {
        coordinate[i] = scale * double(2 * i + 1 - total_point) / double(total_point - 1);
        weight[i] = 0.0;
    }
    for (int i = 0; i < mount; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            weight[i * n + j] += p_weight[j] * dh;
        }
    }
    double max_mag = 0.0;
    if (dimension == 2) {
        for (int j = 0; j < total_point; ++j) {
            for (int i = 0; i < total_point; ++i) {
                Mesh::Cell cell(i + j * total_point);
                cell.position = {coordinate[i], coordinate[j], 0.0};
                cell.volume = weight[i] * weight[j];
                mesh.cells.push_back(cell);
                double cell_mag = cell.position.magnitude();
                if (max_mag < cell_mag) max_mag = cell_mag;
            }
        }
    } else {
        for (int k = 0; k < total_point; ++k) {
            for (int j = 0; j < total_point; ++j) {
                for (int i = 0; i < total_point; ++i) {
                    Mesh::Cell cell(i + (j + k * total_point) * total_point);
                    cell.position = {coordinate[i], coordinate[j], coordinate[k]};
                    cell.volume = weight[i] * weight[j] * weight[k];
                    mesh.cells.push_back(cell);
                    double cell_mag = cell.position.magnitude();
                    if (max_mag < cell_mag) max_mag = cell_mag;
                }
            }
        }
    }
    mesh.update_num();
    mesh.max_cell_magnitude = max_mag;
    return mesh;
}
