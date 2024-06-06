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
void Field<Scalar>::set_zero() {
#pragma omp parallel for shared(len, values, mesh_ptr) default(none)
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = 0.0;
    }
}

template<>
void Field<Vector>::set_zero() {
#pragma omp parallel for shared(len, values, mesh_ptr) default(none)
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = {0.0, 0.0, 0.0};
    }
}

template<>
Field<Vector> Field<Scalar>::gradient(bool _switch) {
    Field<Vector> result(*mesh_ptr, cell_field_flag);
    if (not _switch) return result;
    // openMP - start
#pragma omp parallel for shared(len, values, mesh_ptr, result) default(none)
    for (auto &cell: mesh_ptr->cells) {
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

/// MPI func
void MPI::ReduceAll(Field<MESO::Scalar> &local, Field<MESO::Scalar> &global) {
    for (int i = 0; i < global.len; ++i) {
        MPI_Allreduce(&(local[i]), &(global[i]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

void MPI::ReduceAll(Field<MESO::Vector> &local, Field<MESO::Vector> &global) {
    for (int i = 0; i < global.len; ++i) {
        MPI_Allreduce(&(local[i].x), &(global[i].x), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(local[i].y), &(global[i].y), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(local[i].z), &(global[i].z), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

/// Residual
template<>
Scalar Solver::residual(Field<Scalar> &_old, Field<Scalar> &_new) {
    auto &mesh = *_old.get_mesh();
    ScalarList result_omp(MPI::omp_num, 0.0);
#pragma omp parallel for shared(result_omp, mesh, _old, _new) default(none)
    for (auto &cell: mesh.cells) {
#pragma omp critical
        result_omp[omp_get_team_num()] += std::abs((_new[cell.id] - _old[cell.id]) / _old[cell.id])
                                          * (cell.volume / mesh.total_volume);
    }
    Scalar result = 0.0;
    for (auto it: result_omp) {
        result += it;
    }
    _old = _new;
    return result;
}


template<>
Vector Solver::residual(Field<Vector> &_old, Field<Vector> &_new) {
    auto &mesh = *_old.get_mesh();
    VectorList result_omp(MPI::omp_num, {0.0, 0.0, 0.0});
#pragma omp parallel for shared(result_omp, mesh, _old, _new) default(none)
    for (auto &cell: mesh.cells) {
        auto &new_vec = _new[cell.id];
        auto &old_vec = _old[cell.id];
        Vector vec = {
                std::fabs((new_vec.x - old_vec.x) / old_vec.x),
                std::fabs((new_vec.y - old_vec.y) / old_vec.y),
                std::fabs((new_vec.z - old_vec.z) / old_vec.z)
        };
#pragma omp critical
        result_omp[omp_get_team_num()] += vec * (cell.volume / mesh.total_volume);
    }
    Vector result(0.0, 0.0, 0.0);
    for (auto it: result_omp) {
        result += it;
    }
    _old = _new;
    return result;
}
