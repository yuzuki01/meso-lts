#include "solver/solver.h"


using namespace MESO;


template<>
Field<Scalar>::Field(Mesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    len = (flag == cell_field_flag) ? mesh.NCELL : mesh.NFACE;
    values.resize(len, 0.0);
}

template<>
Field<Vector>::Field(Mesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    len = (flag == cell_field_flag) ? mesh.NCELL : mesh.NFACE;
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
void Field<Scalar>::set_zero() {
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = 0.0;
    }
}

template<>
void Field<Vector>::set_zero() {
    for (int i = 0; i < len; ++i) {
        *(values.begin() + i) = {0.0, 0.0, 0.0};
    }
}

template<>
Field<Vector> Field<Scalar>::gradient(bool _switch) {
    Field<Vector> result(*mesh_ptr, cell_field_flag);
    if (not _switch) return result;
    for (auto &cell: mesh_ptr->cells) {
        Vector Sfr(0.0, 0.0, 0.0);
        for (int j = 0; j < cell.least_square.neighbor_num; ++j) {
            int &neighbor_id = cell.neighbors[j];
            Sfr += (values[neighbor_id] - values[cell.id]) * cell.least_square.weight[j] * cell.least_square.dr[j];
        }
        result[cell.id] = {cell.least_square.Cx * Sfr, cell.least_square.Cy * Sfr, cell.least_square.Cz * Sfr};
    }
    return result;
}

template<>
void Field<Scalar>::output(const std::string &file_name) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    std::stringstream ss;
    ss << file_name << ".np.dat";
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << ss.str() << std::endl;
        fp.close();
        return;
    }
    const int DATA_PRECISION = 18;
    for (auto &it: values) {
        fp << std::setprecision(DATA_PRECISION) << it << std::endl;
    }
    // fp close
    fp.close();
}

template<>
void Field<Vector>::output(const std::string &file_name) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    std::stringstream ss;
    ss << file_name << ".np.dat";
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << ss.str() << std::endl;
        fp.close();
        return;
    }
    const int DATA_PRECISION = 18;
    for (auto &it: values) {
        fp << std::setprecision(DATA_PRECISION) << it.x << " "
           << std::setprecision(DATA_PRECISION) << it.y << " "
           << std::setprecision(DATA_PRECISION) << it.z << std::endl;
    }
    // fp close
    fp.close();
}

/// MPI
void MESO::MPI::AllReduce(Field<MESO::Scalar> &local, Field<MESO::Scalar> &global) {
    MPI_Allreduce(local.values.data(), global.values.data(), global.len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(Field<MESO::Vector> &local, Field<MESO::Vector> &global) {
    MPI_Allreduce(local.values.data(), global.values.data(), global.len, UDF::MPI_Vector, UDF::MPI_VectorSum,
                  MPI_COMM_WORLD);
}

/// Residual
Scalar Solver::residual(Field<Scalar> &_old, Field<Scalar> &_new) {
    auto &mesh = *_old.get_mesh();
    Scalar result = 0.0;
    for (auto &cell: mesh.cells) {
        result += std::abs((_new[cell.id] - _old[cell.id]) / _old[cell.id])
                          * (cell.volume / mesh.total_volume);
    }
    _old = _new;
    return result;
}

Vector Solver::residual(Field<Vector> &_old, Field<Vector> &_new) {
    auto &mesh = *_old.get_mesh();
    Vector result(0.0, 0.0, 0.0);
    for (auto &cell: mesh.cells) {
        auto &new_vec = _new[cell.id];
        auto &old_vec = _old[cell.id];
        Vector vec = {
                std::fabs((new_vec.x - old_vec.x) / old_vec.x),
                std::fabs((new_vec.y - old_vec.y) / old_vec.y),
                std::fabs((new_vec.z - old_vec.z) / old_vec.z)
        };
        result += vec * (cell.volume / mesh.total_volume);
    }
    _old = _new;
    return result;
}
