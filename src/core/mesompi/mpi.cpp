#include "core/core.h"


using namespace MESO;

int MPI::node_num = 1;
int MPI::rank = 0;


void MPI::Initialize(int *p_argc, char ***p_argv) {
    MPI_Init(p_argc, p_argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &node_num);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    std::cout << "Spawn node: " << node_name << " , omp max-num: " << omp_get_max_threads()
              << ", mpi rank: " << rank << " out of " << node_num << " nodes" << std::endl;
}

void MPI::Finalize() {
    MPI_Finalize();
}

MPI::MPI_Task MPI::get_task_distribution(int total_num) {
    int taskPerNode = total_num / node_num;
    int extraTask = total_num % node_num;
    MPI_Task result(MPI::node_num);
    for (int i = 0; i < MPI::node_num; ++i) {
        auto &task = result[i];
        task.start = rank * taskPerNode + (rank < extraTask ? rank : extraTask);
        task.size = taskPerNode + (rank < extraTask ? 1 : 0);
    }
    return result;
}

void MPI::ReduceAndBcast(MESO::Scalar local, MESO::Scalar &global) {
    MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Bcast(&global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPI::ReduceAndBcast(const MESO::Vector &local, MESO::Vector &global) {
    Vector local_sum = local;

    MPI_Reduce(&local_sum.x, &global.x, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum.y, &global.y, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum.z, &global.z, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Bcast(&global.x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global.y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global.z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
