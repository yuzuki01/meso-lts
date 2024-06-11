#include "core/core.h"


using namespace MESO;

int MPI::process_num = 1;
int MPI::rank = 0;

MPI_Datatype MPI::UDF::MPI_Vector;
MPI_Datatype MPI::UDF::type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
int MPI::UDF::block_len[3] = {1, 1, 1};
MPI_Aint MPI::UDF::disp[3] = {offsetof(MESO::Vector, x), offsetof(MESO::Vector, y), offsetof(MESO::Vector, z)};
MPI_Op MPI::UDF::MPI_VectorSum;


void MPI::UDF::vector_sum(void *invec, void *inoutvec, const int *len, MPI_Datatype *datatype) {
    auto *in = (MESO::Vector *) invec;
    auto *inout = (MESO::Vector *) inoutvec;
    for (int i = 0; i < *len; i++) {
        inout[i] += in[i];
    }
}

void MPI::UDF::MPI_UDF_VectorReduce() {
    MPI_Type_create_struct(3, block_len, disp, type, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    MPI_Op_create((MPI_User_function *) vector_sum, 1, &MPI_VectorSum);
}


void MPI::Initialize(int *p_argc, char ***p_argv) {
    MPI_Init(p_argc, p_argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    std::cout << "Spawn node: " << node_name << " , mpi rank: " << rank << " out of " << process_num << " nodes." << std::endl;

    /// 自定义
    UDF::MPI_UDF_VectorReduce();

    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI::Finalize() {
    /// 清理
    MPI_Op_free(&UDF::MPI_VectorSum);
    MPI_Type_free(&UDF::MPI_Vector);

    MPI_Finalize();
}

MPI::MPI_Task MPI::get_task_distribution(int total_num) {
    int taskPerNode = total_num / process_num;
    int extraTask = total_num % process_num;
    MPI_Task result(MPI::process_num);
    for (int i = 0; i < MPI::process_num; ++i) {
        auto &task = result[i];
        task.start = i * taskPerNode + (i < extraTask ? i : extraTask);
        task.size = taskPerNode + (i < extraTask ? 1 : 0);
    }
    return result;
}

void MPI::Bcast(MESO::Scalar &global) {
    MPI_Bcast(&global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPI::Bcast(MESO::Vector &global) {
    MPI_Bcast(&global, 1, UDF::MPI_Vector, 0, MPI_COMM_WORLD);
}

void MPI::AllReduce(MESO::Scalar local, MESO::Scalar &global) {
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MPI::AllReduce(const MESO::Vector &local, MESO::Vector &global) {
    MPI_Allreduce(&local, &global, 1, UDF::MPI_Vector, UDF::MPI_VectorSum, MPI_COMM_WORLD);
}
