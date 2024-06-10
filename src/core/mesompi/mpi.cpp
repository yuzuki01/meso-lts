#include "core/core.h"


using namespace MESO;

int MPI::node_num = 1;
int MPI::rank = 0;
int MPI::omp_num = 1;


void MPI::Initialize(int *p_argc, char ***p_argv, int omp_threads_num) {
    MPI_Init(p_argc, p_argv);

    omp_set_num_threads(omp_threads_num);
    MPI::omp_num = omp_threads_num;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &node_num);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    for (int i = 0; i < node_num; ++i) {
        if (rank == i) {
#pragma omp master
            {
                std::cout << "Spawn node: " << node_name << " , mpi rank: " << rank
                          << " out of " << node_num << " nodes." << std::endl;
            }
#pragma omp parallel shared(node_name, omp_num, std::cout) default(none)
            {
#pragma omp critical
                std::cout << "\tnode: " << node_name << " thread: " << omp_get_thread_num()
                          << " out of " << omp_num << " threads." << std::endl;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
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
        task.start = i * taskPerNode + (i < extraTask ? i : extraTask);
        task.size = taskPerNode + (i < extraTask ? 1 : 0);
    }
    return result;
}

void MPI::Bcast(MESO::Scalar &global) {
    MPI_Bcast(&global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPI::Bcast(MESO::Vector &global) {
    MPI_Bcast(&global.x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global.y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global.z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPI::AllReduce(MESO::Scalar local, MESO::Scalar &global) {
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MPI::AllReduce(const MESO::Vector &local, MESO::Vector &global) {
    MPI_Allreduce(&local.x, &global.x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local.y, &global.y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local.z, &global.z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
