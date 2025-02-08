#include "core/core.h"


namespaceMESO


template<>
void MPI::SendWithStatus(const Label &var, const Label &target, const Label &tag) {
    MPI_Send(&var, 1, MPI_INT, target, tag, MPI_COMM_WORLD);
}

template<>
void MPI::SendWithStatus(const Scalar &var, const Label &target, const Label &tag) {
    MPI_Send(&var, 1, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
}

template<>
void MPI::SendWithStatus(const Vector &var, const Label &target, const Label &tag) {
    MPI_Send(&var, 1, DataType::MPI_VECTOR, target, tag, MPI_COMM_WORLD);
}
