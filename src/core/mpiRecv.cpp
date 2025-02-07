#include "core/core.h"


namespaceMESO


template<>
void MPI::RecvWithStatus(Label &var, const Label &target, Label &tag) {
    MPI_Status status;
    MPI_Recv(&var, 1, MPI_INT, target, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
}

template<>
void MPI::RecvWithStatus(Scalar &var, const Label &target, Label &tag) {
    MPI_Status status;
    MPI_Recv(&var, 1, MPI_DOUBLE, target, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
}

template<>
void MPI::RecvWithStatus(Vector &var, const Label &target, Label &tag) {
    MPI_Status status;
    MPI_Recv(&var, 1, UDF::MPI_Vector, target, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
}
