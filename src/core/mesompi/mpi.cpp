#include "core/core.h"


using namespace MESO;

int MESO::MPI::processor_num = 1;
int MESO::MPI::rank = 0;

MPI_Datatype MESO::MPI::UDF::MPI_Vector;
MPI_Op MESO::MPI::UDF::MPI_VectorSum;

/// local
MPI_Datatype vec_type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
int block_len[3] = {1, 1, 1};
MPI_Aint disp[3] = {offsetof(MESO::Vector, x), offsetof(MESO::Vector, y), offsetof(MESO::Vector, z)};


void MESO::MPI::UDF::vector_sum(void *invec, void *inoutvec, const int *len, MPI_Datatype *datatype) {
    auto *in = (MESO::Vector *) invec;
    auto *inout = (MESO::Vector *) inoutvec;
    for (int i = 0; i < *len; i++) {
        inout[i] += in[i];
    }
}

void MESO::MPI::UDF::MPI_UDF_VectorReduce() {
    MPI_Type_create_struct(3, block_len, disp, vec_type, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    MPI_Op_create((MPI_User_function *) vector_sum, 1, &MPI_VectorSum);
}


void MESO::MPI::Initialize(int *p_argc, char ***p_argv) {
    MPI_Init(p_argc, p_argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processor_num);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    std::cout << "Spawn node: " << node_name << " , mpi rank: " << rank << " out of " << processor_num << " nodes." << std::endl;

    /// 自定义
    UDF::MPI_UDF_VectorReduce();

    MPI_Barrier(MPI_COMM_WORLD);
}

void MESO::MPI::Finalize() {
    /// 清理
    MPI_Op_free(&UDF::MPI_VectorSum);
    MPI_Type_free(&UDF::MPI_Vector);

    MPI_Finalize();
}

MESO::MPI::MPI_Task MESO::MPI::get_task_distribution(int total_num) {
    int taskPerNode = total_num / processor_num;
    int extraTask = total_num % processor_num;
    MPI_Task result(MPI::processor_num);
    for (int i = 0; i < MPI::processor_num; ++i) {
        auto &task = result[i];
        task.start = i * taskPerNode + (i < extraTask ? i : extraTask);
        task.size = taskPerNode + (i < extraTask ? 1 : 0);
    }
    return result;
}

void MESO::MPI::Bcast(bool &global) {
    MPI_Bcast(&global, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
}

void MESO::MPI::Bcast(MESO::Scalar &global) {
    MPI_Bcast(&global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MESO::MPI::Bcast(MESO::Vector &global) {
    MPI_Bcast(&global, 1, UDF::MPI_Vector, 0, MPI_COMM_WORLD);
}

void MESO::MPI::Bcast(List<ObjectId> &global) {
    MPI_Bcast(global.data(), int(global.size()), MPI_INT, 0, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(MESO::Scalar local, MESO::Scalar &global) {
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MESO::MPI::AllReduce(const MESO::Vector &local, MESO::Vector &global) {
    MPI_Allreduce(&local, &global, 1, UDF::MPI_Vector, UDF::MPI_VectorSum, MPI_COMM_WORLD);
}
