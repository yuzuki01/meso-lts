#include "mesh/mesh.h"

using namespace MESO;


Field<Vector> fvmMesh::grad(Field<Scalar> &field) {
    fvmMesh::Mesh *mesh_ptr = field.get_mesh();
    Field<Vector> result(mesh_ptr, field.flag);
    switch (field.flag) {
        case cell_field_flag:
        {
            auto &values = field.values;
            for (auto &cell: mesh_ptr->cells) {
                if (cell.partition_id != MPI::rank) break;
                Vector Sfr(0.0, 0.0, 0.0);
                for (int j = 0; j < cell.least_square.neighbor_num; ++j) {
                    int &neighbor_id = cell.neighbors[j];
                    Sfr += (values[neighbor_id] - values[cell.id]) * cell.least_square.weight[j] *
                           cell.least_square.dr[j];
                }
                result[cell.id] = {cell.least_square.Cx * Sfr, cell.least_square.Cy * Sfr, cell.least_square.Cz * Sfr};
            }
        }
            break;
        case face_field_flag:
        case node_field_flag:
            throw std::invalid_argument("Field<Scalar>.gradient() for flag<face> is not supported.");
        default:
            break;
    }
    return result;
}
