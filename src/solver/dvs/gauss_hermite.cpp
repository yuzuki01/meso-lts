#include "solver/solver.h"


typedef std::unordered_map<int, std::vector<double>> ValueListMap;

ValueListMap gh_weight = {
        {3, {1. / 6., 2. / 3., 1. / 6.}}
};
ValueListMap gh_position = {
        {3, {-sqrt(3), 0.0, sqrt(3)}}
};

using namespace MESO;


Mesh::Mesh Solver::generate_gauss_hermite(int dimension, int _n, double RT) {
    int dvs_num = 1;
    Mesh::Mesh mesh;
    for (int i = 0; i < dimension; ++i) {
        dvs_num *= _n;
    }
    auto &gh_p = gh_position[_n];
    auto &gh_w = gh_weight[_n];
    double max_mag = 0.0;
    if (dimension == 2) {
        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < _n; ++j) {
                Mesh::Cell cell(i + j * _n);
                cell.position = sqrt(RT) * Vector{gh_p[i], gh_p[j], 0.0};
                cell.volume = gh_w[i] * gh_w[j];
                mesh.cells.push_back(cell);
                double cell_mag = cell.position.magnitude();
                if (max_mag < cell_mag) max_mag = cell_mag;
            }
        }
    } else {
        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < _n; ++j) {
                for (int k = 0; k < _n; ++k) {
                    Mesh::Cell cell(i + (j + k * _n) * _n);
                    cell.position = sqrt(RT) * Vector{gh_p[i], gh_p[j], gh_p[k]};
                    cell.volume = gh_w[i] * gh_w[j] * gh_w[k];
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
