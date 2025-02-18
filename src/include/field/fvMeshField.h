#ifndef MESO_FVMESHFIELD_H
#define MESO_FVMESHFIELD_H

namespace MESO::fvm {

    enum {
        LEAST_SQUARE, GREEN_GAUSS
    };

    volVectorField grad(const volScalarField &x, Label defaultMethod=LEAST_SQUARE);

    volVectorField grad(const volScalarField &x, const List<Scalar> &adjData, Label defaultMethod=LEAST_SQUARE);

    template<typename DataType, Label PatchType>
    List<DataType> processorCommAdjData(const BasicField<DataType, PatchType> &x);

    template<typename ValueType>
    List<ValueType> processorCommAllData(const volField<ValueType> &x);
}


template<typename ValueType>
List<ValueType> fvm::processorCommAllData(const volField<ValueType> &x) {
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partitionCellPatch();
    List<ValueType> data(mesh.cellNum());
    forAll(x.index(), ii) {
        data[x.index()[ii]] = x.values()[ii];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 1; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &cell = mesh.cell(patch[pi]);
            if (MPI::rank == cell.rank()) {
                // Send
                MPI::SendWithStatus(data[cell.id()], MPI::mainRank, cell.id());
            } else if (MPI::rank == MPI::mainRank) {
                // Recv
                Label tag = -1;
                ValueType var;
                MPI::RecvWithStatus(var, cell.rank(), tag);
                data[tag] = var;
            }   // Other rank
        }
    }
    MPI::Barrier();
    return data;
}

#endif //MESO_FVMESHFIELD_H
