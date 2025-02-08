#include "core/core.h"


using namespace MESO;


template<>
void MPI::Bcast(Label &global, const Label& root) {
    MPI_Bcast(&global, 1, MPI_INT, root, MPI_COMM_WORLD);
}

template<>
void MPI::Bcast(Scalar &global, const MESO::Label &root) {
    MPI_Bcast(&global, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

template<>
void MPI::Bcast(Vector &global, const MESO::Label &root) {
    MPI_Bcast(&global, 1, DataType::MPI_VECTOR, root, MPI_COMM_WORLD);
}

template<>
void MPI::Bcast(List<Label> &global, const MESO::Label &root) {
    MPI_Bcast(global.data(), Label(global.size()), MPI_INT, root, MPI_COMM_WORLD);
}

template<>
void MPI::Bcast(List<Scalar> &global, const MESO::Label &root) {
    MPI_Bcast(global.data(), Label(global.size()), MPI_DOUBLE, root, MPI_COMM_WORLD);
}

template<>
void MPI::Bcast(List<Vector> &global, const MESO::Label &root) {
    MPI_Bcast(global.data(), Label(global.size()), DataType::MPI_VECTOR, root, MPI_COMM_WORLD);
}
