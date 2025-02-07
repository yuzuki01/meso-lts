#include "field/field.h"


namespaceMESO


List<Scalar> fvm::processorCommAdjData(const volScalarField &x) {
    List<Scalar> data(x.mesh().cellNum());
    const auto& mesh = x.mesh();
    const auto& partition = mesh.partition();
    forAll(x.index(), ci) {
        data[x.index()[ci]] = x[ci];
    }
    if (MPI::processorNum == 1) return data;
    // Send and Recv
    for (int rank = 0; rank < MPI::processorNum; ++rank) {
        const auto& patch = partition[rank];
        forAll(patch.group(), pi) {
            const auto& cell = mesh.cell(patch[pi]);
            forAll(cell.neighbors(), ni) {
                const auto& neighbor = mesh.cell(cell.neighbors()[ni]);
                if (cell.partition() == neighbor.partition()) continue; // Skip
                if (MPI::rank == cell.partition()) {
                    // Send
                    MPI::SendWithStatus(data[cell.id()], neighbor.partition(), neighbor.id());
                } else if (MPI::rank == neighbor.partition()) {
                    // Recv
                    Label tag = -1;
                    MPI::RecvWithStatus(data[neighbor.id()], cell.id(), tag);
                    if (tag != neighbor.id()) {
                        logger.warn << "fvm::processorCommAdjData() receive wrong value" << std::endl;
                    }
                }   // Next
            }
        }
        MPI::Barrier();
    }
    return data;
}


volVectorField gradLeastSquare(const volScalarField &x) {
    auto values = fvm::processorCommAdjData(x);
    const auto &mesh = x.mesh();
    volVectorField grad_(mesh);
    for (const auto &ci : x.index()) {
        const auto &cell = mesh.cell(ci);
        const auto &leastSquare = mesh.leastSquare()[cell.id()];
        Vector Sfr(0.0, 0.0, 0.0);
        forAll(cell.neighbors(), nci) {
            const auto &neighbor = mesh.cell(nci);
            Sfr += (leastSquare.weight[nci]
                    * (values[neighbor.id()] - values[cell.id()]))
                   * leastSquare.dr[nci];
        }
        grad_[ci] = {leastSquare.Cx * Sfr, leastSquare.Cy * Sfr, leastSquare.Cz * Sfr};
    }
    return grad_;
}

volVectorField gradGreenGauss(const volScalarField &x) {
    const auto &mesh = x.mesh();
    volVectorField grad_(mesh);

    return grad_;
}

Label fvm::gradComputationMethod = fvm::LEAST_SQUARE;    // default

/**
 * =========================================================
 * ------------------ FVM Field Functions ------------------
 * =========================================================
 **/


void fvMesh::initLeastSquare() {
    if (not leastSquare_.empty()) return;
    leastSquare_.reserve(NCELL);
    for (const auto &cell: cells_) {        // after fvMesh::partition()
        const auto neiSize = static_cast<Label>(cell.neighbors().size());
        LeastSquare ls(neiSize);
        Scalar Sxx, Sxy, Sxz, Syy, Syz, Szz;
        Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
        for (int i = 0; i < neiSize; i++) {
            auto &neiCell = cells_[cell.neighbors()[i]];
            Vector dr = neiCell.C() - cell.C();
            Scalar wi = 1.0 / (dr * dr);
            ls.dr[i] = dr;
            ls.weight[i] = wi;
            // sum
            Sxx += wi * (dr.x * dr.x);
            Sxy += wi * (dr.x * dr.y);
            Sxz += wi * (dr.x * dr.z);
            Syy += wi * (dr.y * dr.y);
            Syz += wi * (dr.y * dr.z);
            Szz += wi * (dr.z * dr.z);
        }
        if (dimension() == 2) {
            Scalar FM = Sxx * Syy - Sxy * Sxy;
            ls.Cx = {Syy / FM, -Sxy / FM, 0.0};
            ls.Cy = {-Sxy / FM, Sxx / FM, 0.0};
            ls.Cz = {0.0, 0.0, 0.0};
        } else {
            Scalar FM = -Sxz * Sxz * Syy + 2.0 * Sxy * Sxz * Syz - Sxx * Syz * Syz - Sxy * Sxy * Szz +
                        Sxx * Syy * Szz;
            Scalar SxyyzNxzyy = Sxy * Syz - Sxz * Syy;
            Scalar SxzyzNxyzz = Sxz * Syz - Sxy * Szz;
            Scalar SyyzzNyzyz = Syy * Szz - Syz * Syz;
            Scalar SxxzzNxzxz = Sxx * Szz - Sxz * Sxz;
            Scalar SxyxzNxxyz = Sxy * Sxz - Sxx * Syz;
            Scalar SxxyyNxyxy = Sxx * Syy - Sxy * Sxy;
            ls.Cx = {
                    SyyzzNyzyz / FM, SxzyzNxyzz / FM, SxyyzNxzyy / FM
            };
            ls.Cy = {
                    SxzyzNxyzz / FM, SxxzzNxzxz / FM, SxyxzNxxyz / FM
            };
            ls.Cz = {
                    SxyyzNxzyy / FM, SxyxzNxxyz / FM, SxxyyNxyxy / FM
            };
        }
        leastSquare_.push_back(ls);
    }
}

const List<LeastSquare> &fvMesh::leastSquare() const {
    return leastSquare_;
}

volVectorField fvm::grad(const volScalarField &x) {
    switch (fvm::gradComputationMethod) {
        case fvm::LEAST_SQUARE:
            return gradLeastSquare(x);
        case fvm::GREEN_GAUSS:
            return gradGreenGauss(x);
        default:
            logger.error << "fvm::grad() caught unexpected value" << std::endl;
            FATAL_ERROR_THROW;
    }
}
