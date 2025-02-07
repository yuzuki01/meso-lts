#include "mesh/mesh.h"

using namespace MESO;

Vector fvmMesh::grad(Field<MESO::Scalar> &field, MESO::ObjectId element_id) {
    fvmMesh::Mesh *mesh_ptr = field.get_mesh();
    switch (field.flag) {
        case cell_field_flag: {
            auto &cell = mesh_ptr->cells[element_id];
            if (cell.partition_id != MPI::rank) return {0.0, 0.0, 0.0};
            auto &values = field.values;
            Vector Sfr(0.0, 0.0, 0.0);
            for (int j = 0; j < cell.least_square.neighbor_num; ++j) {
                int &neighbor_id = cell.neighbors[j];
                Sfr += (values[neighbor_id] - values[cell.id]) * cell.least_square.dr[j];
            }
            return {cell.least_square.Cx * Sfr, cell.least_square.Cy * Sfr, cell.least_square.Cz * Sfr};
        }
        default:
            throw std::invalid_argument("Field<Scalar>.gradient() for flag<face> is not supported.");
    }
}

Field<Vector> fvmMesh::grad(Field<Scalar> &field) {
    fvmMesh::Mesh *mesh_ptr = field.get_mesh();
    Field<Vector> result(mesh_ptr, field.flag);
    switch (field.flag) {
        case cell_field_flag: {
            auto &values = field.values;
            for (auto &cell: mesh_ptr->cells) {
                if (cell.partition_id != MPI::rank) break;
                Vector Sfr(0.0, 0.0, 0.0);
                for (int j = 0; j < cell.least_square.neighbor_num; ++j) {
                    int &neighbor_id = cell.neighbors[j];
                    auto wj = cell.least_square.weight[j];
                    Sfr += wj * (values[neighbor_id] - values[cell.id]) * cell.least_square.dr[j];
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

Scalar fvmMesh::interp_IDW(Field<MESO::Scalar> &field, MESO::ObjectId element_id) {
    fvmMesh::Mesh *mesh_ptr = field.get_mesh();
    switch (field.flag) {
        case cell_field_flag: {
            auto &cell = mesh_ptr->cells[element_id];
            if (cell.partition_id != MPI::rank) return field[cell.id];
            auto &values = field.values;
            Scalar sum_wf = 0.0, sum_w = 0.0;
            for (int i = 0; i < cell.least_square.neighbor_num; ++i) {
                auto neighbor = cell.neighbors[i];
                auto &dr = cell.least_square.dr[i];
                auto wi = 1.0 / (dr * dr);
                sum_w += wi;
                sum_wf += wi * values[neighbor];
            }
            return sum_wf / sum_w;
        }
        default:
            throw std::invalid_argument("Field<Scalar>.gradient() for flag<face> is not supported.");
    }
}

Vector fvmMesh::interp_IDW(Field<MESO::Vector> &field, MESO::ObjectId element_id) {
    fvmMesh::Mesh *mesh_ptr = field.get_mesh();
    switch (field.flag) {
        case cell_field_flag: {
            auto &cell = mesh_ptr->cells[element_id];
            if (cell.partition_id != MPI::rank) return field[cell.id];
            auto &values = field.values;
            Vector sum_wf(0.0, 0.0, 0.0);
            Scalar sum_w = 0.0;
            for (int i = 0; i < cell.least_square.neighbor_num; ++i) {
                auto neighbor = cell.neighbors[i];
                auto &dr = cell.least_square.dr[i];
                auto wi = 1.0 / (dr * dr);
                sum_w += wi;
                sum_wf += wi * values[neighbor];
            }
            return sum_wf / sum_w;
        }
        default:
            throw std::invalid_argument("Field<Scalar>.gradient() for flag<face> is not supported.");
    }
}
