#include "field/field.h"


namespaceMESO


volVectorField gradLeastSquare(const volScalarField &x) {
    auto values = fvm::processorCommAdjData(x);
    const auto &mesh = x.mesh();
    volVectorField grad_(mesh);
    forConstRef(ci, x.index()) {
        const auto &cell = mesh.cell(ci);
        const auto &leastSquare = mesh.leastSquare()[cell.id()];
        Vector Sfr(0.0, 0.0, 0.0);
        forAll(cell.neighbors(), ni) {
            const auto &neighbor = mesh.cell(cell.neighbors()[ni]);
            Sfr += (leastSquare.weight[ni]
                    * (values[cell.id()] - values[neighbor.id()]))
                   * leastSquare.dr[ni];
        }
        grad_.values()[cell.idOnPartition()] = {leastSquare.Cx * Sfr, leastSquare.Cy * Sfr, leastSquare.Cz * Sfr};
    }
    return grad_;
}


volVectorField gradGreenGauss(const volScalarField &x, const Label iterTimes) {
    auto valuesCell = fvm::processorCommAdjData(x);
    const auto &mesh = x.mesh();
    volVectorField grad_(mesh);
    surfScalarField x_f(x.mesh());
    // Step.1 - predict x_f*
    forConstRef(fi, x_f.index()) {
        const auto &face = mesh.face(fi);
        const auto &owner = mesh.cell(face.owner());
        const auto &neighbor = mesh.cell(face.neighbor());
        auto drMag_of = (face.C() - owner.C()).magnitude();
        auto drMag_nf = (face.C() - neighbor.C()).magnitude();
        Scalar w_o = drMag_of / (drMag_of + drMag_nf);
        x_f.values()[face.idOnPartition()] = w_o * valuesCell[owner.id()] + (1 - w_o) * valuesCell[neighbor.id()];
    }
    // Step.2 - compute grad_c
    auto valuesFace = fvm::processorCommAdjData(x_f);
    forConstRef(ci, x.index()) {
        const auto &cell = mesh.cell(ci);
        Vector Sf(0, 0, 0);
        forConstRef(fi, cell.faces()) {
            const auto &face = mesh.face(fi);
            auto Snv = (face.owner() == cell.id()) ? (face.Snv()) : (-face.Snv());
            Sf += valuesFace[face.id()] * Snv;
        }
        grad_.values()[cell.idOnPartition()] = Sf / cell.V();
    }
    List<Vector> valuesGrad;
    for (int i = 0; i < iterTimes; ++i) {
        // Step.3 - compute x_f
        valuesGrad = fvm::processorCommAdjData(grad_);
        forConstRef(fi, x_f.index()) {
            const auto &face = mesh.face(fi);
            const auto &owner = mesh.cell(face.owner());
            const auto &neighbor = mesh.cell(face.neighbor());
            auto dr_of = face.C() - owner.C();
            auto dr_nf = face.C() - neighbor.C();
            auto drMag_of = dr_of.magnitude();
            auto drMag_nf = dr_nf.magnitude();
            Scalar w_o = drMag_of / (drMag_of + drMag_nf);
            x_f.values()[face.idOnPartition()] = w_o * (valuesCell[owner.id()] + valuesGrad[owner.id()] * dr_of)
                                                 + (1 - w_o) * (valuesCell[neighbor.id()] + valuesGrad[neighbor.id()] * dr_nf);
        }
        // Step.4 - compute grad_c again
        valuesFace = fvm::processorCommAdjData(x_f);
        forConstRef(ci, x.index()) {
            const auto &cell = mesh.cell(ci);
            Vector Sf(0, 0, 0);
            forConstRef(fi, cell.faces()) {
                const auto &face = mesh.face(fi);
                auto Snv = (face.owner() == cell.id()) ? (face.Snv()) : (-face.Snv());
                Sf += valuesFace[face.id()] * Snv;
            }
            grad_.values()[cell.idOnPartition()] = Sf / cell.V();
        }
    }
    return grad_;
}


/**
 * =========================================================
 * ------------------ FVM Field Functions ------------------
 * =========================================================
 **/


void fvMesh::initLeastSquare() {
    leastSquare_.clear();
    leastSquare_.reserve(NCELL);
    forConstRef(cell, cells_) {     // after fvMesh::partition()
        const auto neighborNum = static_cast<Label>(cell.neighbors().size());
        LeastSquare ls(neighborNum);
        Scalar Sxx, Sxy, Sxz, Syy, Syz, Szz;
        Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
        forAll(cell.neighbors(), i) {
            const auto &neighbor = cells_[cell.neighbors()[i]];
            Vector dr = cell.C() - neighbor.C();
            Scalar wi = 1.0 / magSqr(dr);
            ls.dr[i] = dr;
            ls.weight[i] = wi;
            // sum
            Sxx += (dr.x * dr.x) * wi;
            Sxy += (dr.x * dr.y) * wi;
            Sxz += (dr.x * dr.z) * wi;
            Syy += (dr.y * dr.y) * wi;
            Syz += (dr.y * dr.z) * wi;
            Szz += (dr.z * dr.z) * wi;
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

volVectorField fvm::grad(const volScalarField &x, Label defaultMethod) {
    switch (defaultMethod) {
        case fvm::LEAST_SQUARE:
            return gradLeastSquare(x);
        case fvm::GREEN_GAUSS:
            return gradGreenGauss(x, 3);
        default:
            logger.error << "fvm::grad() caught unexpected value" << std::endl;
            FATAL_ERROR_THROW;
    }
}
