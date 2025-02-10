#include "field/field.h"


namespaceMESO


template<>
List<Scalar> fvm::processorCommAdjData(const volScalarField &x) {
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partitionCellPatch();
    List<Scalar> data(mesh.cellNum());
    forAll(x.index(), ii) {
        data[x.index()[ii]] = x.values()[ii];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &cell = mesh.cell(patch[pi]);
            forConstRef(ni, cell.neighbors()) {
                const auto &neighbor = mesh.cell(ni);
                if (cell.rank() == neighbor.rank()) continue;
                if (MPI::rank == cell.rank()) {
                    // Send
                    MPI::SendWithStatus(data[cell.id()], neighbor.rank(), cell.id());
                } else if (MPI::rank == neighbor.rank()) {
                    // Recv
                    Label tag = -1;
                    Scalar var;
                    MPI::RecvWithStatus(var, cell.rank(), tag);
                    data[tag] = var;
                }   // Other rank
            }
        }
    }
    MPI::Barrier();
    return data;
}


template<>
List<Scalar> fvm::processorCommAdjData(const surfScalarField &x) {
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partitionFacePatch();
    List<Scalar> data(mesh.faceNum());
    forAll(x.index(), ii) {
        data[x.index()[ii]] = x.values()[ii];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &face = mesh.face(patch[pi]);
            const auto &owner = mesh.cell(face.owner());
            const auto &neighbor = mesh.cell(face.neighbor());
            if (face.rank() == neighbor.rank()) continue;
            if (MPI::rank == owner.rank()) {
                // Send
                MPI::SendWithStatus(data[face.id()], neighbor.rank(), face.id());
            } else if (MPI::rank == neighbor.rank()) {
                // Recv
                Label tag = -1;
                Scalar var;
                MPI::RecvWithStatus(var, face.rank(), tag);
                data[tag] = var;
            }   // Other rank
        }
    }
    MPI::Barrier();
    return data;
}


template<>
List<Vector> fvm::processorCommAdjData(const volVectorField &x) {
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partitionCellPatch();
    List<Vector> data(mesh.cellNum());
    forAll(x.index(), ii) {
        data[x.index()[ii]] = x.values()[ii];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &cell = mesh.cell(patch[pi]);
            forConstRef(ni, cell.neighbors()) {
                const auto &neighbor = mesh.cell(ni);
                if (cell.rank() == neighbor.rank()) continue;
                if (MPI::rank == cell.rank()) {
                    // Send
                    MPI::SendWithStatus(data[cell.id()], neighbor.rank(), cell.id());
                } else if (MPI::rank == neighbor.rank()) {
                    // Recv
                    Label tag = -1;
                    Vector var;
                    MPI::RecvWithStatus(var, cell.rank(), tag);
                    data[tag] = var;
                }   // Other rank
            }
        }
    }
    MPI::Barrier();
    return data;
}


template<>
List<Vector> fvm::processorCommAdjData(const surfVectorField &x) {
    const auto &mesh = x.mesh();
    const auto &partition = mesh.partitionFacePatch();
    List<Vector> data(mesh.faceNum());
    forAll(x.index(), ii) {
        data[x.index()[ii]] = x.values()[ii];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto &patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto &face = mesh.face(patch[pi]);
            const auto &owner = mesh.cell(face.owner());
            const auto &neighbor = mesh.cell(face.neighbor());
            if (face.rank() == neighbor.rank()) continue;
            if (MPI::rank == owner.rank()) {
                // Send
                MPI::SendWithStatus(data[face.id()], neighbor.rank(), face.id());
            } else if (MPI::rank == neighbor.rank()) {
                // Recv
                Label tag = -1;
                Vector var;
                MPI::RecvWithStatus(var, face.rank(), tag);
                data[tag] = var;
            }   // Other rank
        }
    }
    MPI::Barrier();
    return data;
}
