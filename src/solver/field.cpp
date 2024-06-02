#include "solver/solver.h"


using namespace MESO;


template<>
Field<Scalar>::Field(Mesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    len = (flag == cell_field_flag) ? mesh.NCELL : int(mesh.faces.size());
    values.resize(len, 0.0);
}

template<>
Field<Vector>::Field(Mesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    len = (flag == cell_field_flag) ? mesh.NCELL : int(mesh.faces.size());
    values.resize(len, Vector(0.0, 0.0, 0.0));
}

template<>
Field<Scalar>::Field(const Field<Scalar> &other) : mesh_ptr(other.get_mesh()) {
    flag = other.flag;
    values = other.values;
}

template<>
Field<Vector>::Field(const Field<Vector> &other) : mesh_ptr(other.get_mesh()) {
    flag = other.flag;
    values = other.values;
}

template<>
Field<Scalar> Field<Vector>::heft(int id) {
    Field<Scalar> result(*mesh_ptr, flag);
#pragma omp parallel for shared(len, result, id) default(none)
    for (int i = 0; i < len; ++i) {
        result.values[i] = values[i][id];
    }
    return result;
}

template<>
Scalar &Field<Scalar>::operator[](int index) {
    return values[index];
}

template<>
Vector &Field<Vector>::operator[](int index) {
    return values[index];
}

Field<Scalar> Mesh::Mesh::zero_scalar_field(int flag) {
    return {*this, flag};
}

Field<Vector> Mesh::Mesh::zero_vector_field(int flag) {
    return {*this, flag};
}

template<>
void Field<Scalar>::MeshCellValueToField(const std::function<Scalar(Mesh::Cell &)> &func) {
#pragma omp parallel for shared(len, values, mesh_ptr, func) default(none)
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = func(*(mesh_ptr->cells.begin() + i));
    }
}

template<>
void Field<Vector>::MeshCellValueToField(const std::function<Vector(Mesh::Cell &)> &func) {
#pragma omp parallel for shared(len, values, mesh_ptr, func) default(none)
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = func(*(mesh_ptr->cells.begin() + i));
    }
}

template<>
Field<Vector> Field<Scalar>::gradient(bool _switch) {
    Field<Vector> result(*mesh_ptr, cell_field_flag);
    if (not _switch) return result;
    // openMP - start
#pragma omp parallel for shared(len, values, mesh_ptr, result) default(none)
    for (auto &cell : mesh_ptr->cells) {
        Vector Sfr(0.0, 0.0, 0.0);
        for (int j = 0; j < cell.least_square.neighbor_num; ++j) {
            int &neighbor_id = cell.neighbors[j];
            Sfr += (values[neighbor_id] - values[cell.id]) * cell.least_square.weight[j] * cell.least_square.dr[j];
        }
        result[cell.id] = {cell.least_square.Cx * Sfr, cell.least_square.Cy * Sfr, cell.least_square.Cz * Sfr};
    }
    // openMP - end
    return result;
}
