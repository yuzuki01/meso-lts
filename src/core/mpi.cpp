#include "core/core.h"


using namespace MESO;

Label MPI::processorNum = 1;
Label MPI::rank = 0;
MPI_Datatype MPI::UDF::MPI_Vector;

void MPI::Initialize(Label *p_argc, char ***p_argv) {
    MPI_Init(p_argc, p_argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processorNum);

    char nodeName[MPI_MAX_PROCESSOR_NAME];
    Label name_len;
    MPI_Get_processor_name(nodeName, &name_len);

    std::cout << "Spawn node: " << nodeName
              << " , mpi rank: " << rank
              << " out of " << processorNum
              << " nodes." << std::endl;

    /// 自定义
    // functions

    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI::Finalize() {
    /// 清理
    // MPI_Type_free(&UDF::MPI_Vector);

    MPI_Finalize();
}

void MPI::Barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}
