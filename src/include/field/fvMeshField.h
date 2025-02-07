#ifndef MESO_FVMESHFIELD_H
#define MESO_FVMESHFIELD_H

namespace MESO::fvm {

    enum {
        LEAST_SQUARE, GREEN_GAUSS
    };
    extern Label gradComputationMethod;

    volVectorField grad(const volScalarField &x);

    List<Scalar> processorCommAdjData(const volScalarField &x);

    template<typename ValueType>
    List<ValueType> processorCommAllData(const volField<ValueType> &x);
}


template<typename ValueType>
List<ValueType> fvm::processorCommAllData(const volField<ValueType> &x) {
    List<ValueType> data(x.mesh().cellNum());
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partition();
    forAll(x.index(), ci) {
        data[x.index()[ci]] = x[ci];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &cell = mesh.cell(patch[pi]);
            if (MPI::rank == MPI::mainRank) {
                // Recv
                Label tag = -1;
                ValueType var;
                MPI::RecvWithStatus(var, cell.partition(), tag);
                data[tag] = var;
            } else if (MPI::rank == cell.partition()) {
                // Send
                MPI::SendWithStatus(data[cell.id()], MPI::mainRank, cell.id());
            }  // Next
        }
        MPI::Barrier();
    }
    return data;
}

#endif //MESO_FVMESHFIELD_H
