#include "mesh.h"

/// local vars
std::unordered_map<int, std::vector<double>> gauss_hermit_points{
        {3, {-sqrt(3.0), 0.0, sqrt(3.0)}}
};


std::unordered_map<int, std::vector<double>> gauss_hermit_weights{
        {3, {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}}
};

/*************************************************
 * gauss_hermit_points                           *
 *     - GH 型积分的高斯点，乘于 sqrt(RT) 可以得到速度 *
 * gauss_hermit_weights                          *
 *     - GH 型积分的权重                           *
 *************************************************/

MESH::StaticMesh GENERATOR::gauss_hermit(int gauss_point, int dimension, double RT) {
    std::stringstream mesh_name;
    int Q = 1;
    for (int i = 0; i < dimension; i++) Q *= gauss_point;
    mesh_name << "D" << dimension << "Q" << Q;
    MESH::StaticMesh mesh;
    mesh.name = mesh_name.str();
    auto &points = gauss_hermit_points[gauss_point];
    auto &weights = gauss_hermit_weights[gauss_point];
    if (dimension == 2) {
        for (int j = 0; j < gauss_point; j++)
            for (int i = 0; i < gauss_point; i++) {
                /// cell_key from 0 - Q-1
                int cell_key = i + j * gauss_point;
                int inv_cell_key = (gauss_point - i - 1) + (gauss_point - j - 1) * gauss_point;
                mesh.CELLS.emplace_back(
                        sqrt(RT) * Vec3D(points[i], points[j], 0.0),
                        weights[i] * weights[j],
                        cell_key, inv_cell_key, -1, -1
                );
                auto &cell = mesh.get_cell(cell_key);
                mesh.update_max_discrete_velocity(cell.position.magnitude());
            }
    } else if (dimension == 3) {
        for (int k = 0; k < gauss_point; k++)
            for (int j = 0; j < gauss_point; j++)
                for (int i = 0; i < gauss_point; i++) {
                    /// cell_key from 0 - Q-1
                    int cell_key = i + j * gauss_point + k * gauss_point * gauss_point;
                    int inv_cell_key =
                            (gauss_point - i - 1) + (gauss_point - j - 1) * gauss_point +
                            (gauss_point - k - 1) * gauss_point;
                    mesh.CELLS.emplace_back(
                            sqrt(RT) * Vec3D(points[i], points[j], points[k]),
                            weights[i] * weights[j] * weights[k],
                            cell_key, inv_cell_key, -1, -1
                    );
                    auto &cell = mesh.get_cell(cell_key);
                    mesh.update_max_discrete_velocity(cell.position.magnitude());
                }
    }
    mesh.NDIME = dimension;
    mesh.NELEM = int(mesh.CELLS.size());
    mesh.shrink_to_fit();
    std::stringstream ss;
    ss << "Generate DVS mesh <" << mesh_name.str() << ">";
    note_println(ss.str());
    return mesh;
}
