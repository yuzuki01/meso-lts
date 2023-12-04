#include "mesh.h"

/// local vars

std::unordered_map<int, std::vector<double>> newton_cotes_weight = {
        {1, {1./2., 1./2.}},
        {2, {1./6., 4./6., 1./6.}},
        {3, {1./8., 3./8., 3./8., 1./8.}},
        {4, {7./90., 16./45., 2./15., 16./45., 7./90.}},
        {5, {19./288., 25./96., 25./144., 25./144., 25./96., 19./288.}},
        {6, {41./840., 9./35., 9./280., 34./105., 9./280., 9./35., 41./840.}}
};

/*************************************************
 * Newton Cotes                                  *
 *   weight : {n, {a0, a1, .., an}}              *
 *   Integral{f} = Sum{ai * f(xi)} + O(x^(2n-1)) *
 *************************************************/

MESH::StaticMesh GENERATOR::newton_cotes(int n, int mount, int dimension, double scale, double RT) {
    std::stringstream mesh_name;
    int total_point = n * mount + 1;
    auto &point_weight = newton_cotes_weight.at(n);
    mesh_name << "newton-cotes(" << n << ")-" << total_point << "x" << total_point;
    if (dimension == 3) mesh_name << "x" << total_point;
    MESH::StaticMesh mesh;
    mesh.name = mesh_name.str();
    scale = scale * sqrt(2.0 * RT);
    double dh = 2.0 * scale / mount;
    std::vector<double> coordinate(total_point), weight(total_point);
    for (int i = 0; i < total_point; i++) {
        coordinate[i] = scale * double(2 * i + 1 - total_point) / double (total_point - 1);
        weight[i] = 0.0;
    }
    for (int i = 0; i < mount; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            weight[i * n + j] += point_weight[j] * dh;
        }
    }
    if (dimension == 2) {
        for (int j = 0; j < total_point; ++j)
            for (int i = 0; i < total_point; ++i) {
                int cell_key = i + j * total_point;
                int inv_cell_key = (total_point - i - 1) + (total_point - j - 1) * total_point;
                mesh.CELLS.emplace_back(
                        Vec3D(coordinate[i], coordinate[j], 0.0),
                        weight[i] * weight[j],
                        cell_key, inv_cell_key, -1, -1
                        );
                auto &cell = mesh.get_cell(cell_key);
                mesh.update_max_discrete_velocity(cell.position.magnitude());
            }
    } else if (dimension == 3) {
        for (int k = 0; k < total_point; k++)
            for (int j = 0; j < total_point; j++)
                for (int i = 0; i < total_point; i++) {
                    int cell_key = i + j * total_point + k * total_point * total_point;
                    int inv_cell_key =
                            (total_point - i - 1) + (total_point - j - 1) * total_point +
                            (total_point - k - 1) * total_point * total_point;
                    mesh.CELLS.emplace_back(
                            Vec3D(coordinate[i], coordinate[j], coordinate[k]),
                            weight[i] * weight[j] * weight[k],
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
