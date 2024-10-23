#include "solver/solver.h"


using namespace MESO;


template<>
Field<Scalar>::Field(fvmMesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    switch (flag) {
        case cell_field_flag:
            len = mesh.NCELL;
            break;
        case face_field_flag:
            len = mesh.NFACE;
            break;
        case node_field_flag:
            len = mesh.NNODE;
            break;
        default:
            logger.error << "Field<Scalar> caught wrong flag: " << flag << std::endl;
            throw std::invalid_argument("Field<Scalar> caught wrong flag");
    }
    values.resize(len, 0.0);
}

template<>
Field<Scalar>::Field(fvmMesh::Mesh *mesh_ptr, int flag) : mesh_ptr(mesh_ptr), flag(flag) {
    auto &mesh = *mesh_ptr;
    switch (flag) {
        case cell_field_flag:
            len = mesh.NCELL;
            break;
        case face_field_flag:
            len = mesh.NFACE;
            break;
        case node_field_flag:
            len = mesh.NNODE;
            break;
        default:
            logger.error << "Field<Scalar> caught wrong flag: " << flag << std::endl;
            throw std::invalid_argument("Field<Scalar> caught wrong flag");
    }
    values.resize(len, 0.0);
}

template<>
Field<Vector>::Field(fvmMesh::Mesh &mesh, int flag) : mesh_ptr(&mesh), flag(flag) {
    switch (flag) {
        case cell_field_flag:
            len = mesh.NCELL;
            break;
        case face_field_flag:
            len = mesh.NFACE;
            break;
        case node_field_flag:
            len = mesh.NNODE;
            break;
        default:
            logger.error << "Field<Scalar> caught wrong flag: " << flag << std::endl;
            throw std::invalid_argument("Field<Vector> caught wrong flag");
    }
    values.resize(len, Vector(0.0, 0.0, 0.0));
}

template<>
Field<Vector>::Field(fvmMesh::Mesh *mesh_ptr, int flag) : mesh_ptr(mesh_ptr), flag(flag) {
    auto &mesh = *mesh_ptr;
    switch (flag) {
        case cell_field_flag:
            len = mesh.NCELL;
            break;
        case face_field_flag:
            len = mesh.NFACE;
            break;
        case node_field_flag:
            len = mesh.NNODE;
            break;
        default:
            logger.error << "Field<Scalar> caught wrong flag: " << flag << std::endl;
            throw std::invalid_argument("Field<Vector> caught wrong flag");
    }
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

Field<Scalar> fvmMesh::Mesh::zero_scalar_field(int flag) {
    return {*this, flag};
}

Field<Vector> fvmMesh::Mesh::zero_vector_field(int flag) {
    return {*this, flag};
}

template<>
void Field<Scalar>::set_zero() {
    for (int i = 0; i < len; ++i) {
        values[i] = 0.0;
    }
}

template<>
void Field<Vector>::set_zero() {
    for (int i = 0; i < len; ++i) {
        values[i] = {0.0, 0.0, 0.0};
    }
}

template<>
void Field<Scalar>::output(const std::string &file_name) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    fp.open(file_name, std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << file_name << std::endl;
        fp.close();
        return;
    }
    const int DATA_PRECISION = 18;
    for (auto &it: values) {
        fp << std::setprecision(DATA_PRECISION) << it << std::endl;
    }
    fp.close();
}

template<>
void Field<Vector>::output(const std::string &file_name) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    fp.open(file_name, std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << file_name << std::endl;
        fp.close();
        return;
    }
    const int DATA_PRECISION = 18;
    for (auto &it: values) {
        fp << std::setprecision(DATA_PRECISION) << it.x << " "
           << std::setprecision(DATA_PRECISION) << it.y << " "
           << std::setprecision(DATA_PRECISION) << it.z << std::endl;
    }
    fp.close();
}

/// MPI
void MESO::MPI::AllReduce(MESO::List<Scalar> &local, MESO::List<Scalar> &global) {
    if (local.size() != global.size()) global.resize(local.size());
    MPI_Allreduce(local.data(), global.data(), int(global.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(MESO::List<Vector> &local, MESO::List<Vector> &global) {
    if (local.size() != global.size()) global.resize(local.size());
    MPI_Allreduce(local.data(), global.data(), int(global.size()), UDF::MPI_Vector, UDF::MPI_VectorSum, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(Field<MESO::Scalar> &local, Field<MESO::Scalar> &global) {
    MPI_Allreduce(local.values.data(), global.values.data(), global.len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(Field<MESO::Vector> &local, Field<MESO::Vector> &global) {
    MPI_Allreduce(local.values.data(), global.values.data(), global.len, UDF::MPI_Vector, UDF::MPI_VectorSum,
                  MPI_COMM_WORLD);
}

/// Gather
void
MESO::MPI::GatherFieldList(std::vector<Field<Scalar>> &local, MESO::List<Scalar> &global, MESO::ObjectId column_id) {
    int local_size = static_cast<int>(local.size());
    List<ObjectId> recv_counts(MPI::processor_num);
    MPI_Gather(&local_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    List<ObjectId> displacements(MPI::processor_num, 0);
    if (MPI::rank == 0) {
        for (int i = 0; i < MPI::processor_num; ++i) {
            displacements[i] = displacements[i - 1] + recv_counts[i - 1];
        }
    }
    MESO::List<Scalar> gathered_data(displacements[MPI::processor_num - 1] + recv_counts[MPI::processor_num - 1]);
    MESO::List<Scalar> local_list(local_size);
    for (int i = 0; i < local_size; ++i) {
        local_list[i] = local[i][column_id];
    }
    MPI_Gatherv(local_list.data(), local_size, MPI_DOUBLE, gathered_data.data(), recv_counts.data(),
                displacements.data(),
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < static_cast<int>(gathered_data.size()); ++i) {
        global[i] = gathered_data[i];
    }
    MPI_Bcast(global.data(), static_cast<int>(global.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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


/// numpy API
template<>
void MESO::read_np_data<MESO::Scalar>(const std::string &file, Field<MESO::Scalar> &field) {
    MESO::FileReader::BasicReader reader(file);
    auto lines = reader.read_lines();
    auto line_num = static_cast<int>(lines.size());
    for (int i = 0; i < line_num; ++i) {
        auto data = lines[i];
        field[i] = stod(data);
    }
}

template<>
void MESO::read_np_data<MESO::Vector>(const std::string &file, Field<MESO::Vector> &field) {
    MESO::FileReader::BasicReader reader(file);
    auto lines = reader.read_lines();
    auto line_num = static_cast<int>(lines.size());
    for (int i = 0; i < line_num; ++i) {
        auto data = MESO::Utils::split(lines[i]);
        field[i] = {stod(data[0]), stod(data[1]), stod(data[2])};
    }
}
