#include "core/core.h"


using namespace MESO;

int MPI::processor_num = 1;
int MPI::rank = 0;

void MPI::Initialize(int *p_argc, char ***p_argv) {
    MPI_Init(p_argc, p_argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processor_num);

    char node_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(node_name, &name_len);

    std::cout << "Spawn node: " << node_name << " , omp max-num: " << omp_get_max_threads()
              << ", mpi rank: " << rank << " out of " << processor_num << " nodes" << std::endl;
}

void MPI::Finalize() {
    MPI_Finalize();
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
