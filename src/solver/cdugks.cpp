#include "solver/solver.h"

using namespace MESO;
using namespace Solver;


CDUGKS::CDUGKS(const MESO::String &filePath)
        : BasicSolver(filePath),
          DV_(config_.get<String>("DVMProperties",
                                  "mesh",
                                  "constant/dvs.neu"),
              time_,
              true),
          R(config_.get("DVMProperties", "gasConstant", 0.5)),
          TRef(config_.get("DVMProperties", "TemperatureRef", 1.0)),
          muRef(config_.get("DVMProperties", "ViscosityRef", 1.0e-5)),
          rhoVol(mesh_), rhoVolOld(mesh_),
          UVol(mesh_), UVolOld(mesh_),
          rhoSurf(mesh_),
          USurf(mesh_),
          gVol_(DV_.cellNum(), volScalarField(mesh_)),
          gSurf_(DV_.cellNum(), surfScalarField(mesh_)),
          fluxGVol_(DV_.cellNum(), volScalarField(mesh_)) {
    mesh_.info();
}

/// CDUGKS
Scalar CDUGKS::gMaxwell(const Scalar &rho_, const Vector &U_, const Cell &dv_) {
    auto RT2 = 2.0 * R * TRef;
    auto c = dv_.C() - U_;
    return rho_ / pow((M_PI * RT2), mesh_.dimension() * 0.5) * exp(-(c * c) / RT2);
}

template<Label PatchType>
BasicField<Scalar, PatchType> CDUGKS::gMaxwell(const BasicField<Scalar, PatchType> &rho_,
                                               const BasicField<Vector, PatchType> &U_,
                                               const Mesh::Cell &dv_) {
    auto RT2 = 2.0 * R * TRef;
    auto c = dv_.C() - U_;
    return rho_ / pow((M_PI * RT2), mesh_.dimension() * 0.5) * exp(-(c * c) / RT2);
}

template<Label PatchType>
BasicField<Scalar, PatchType> CDUGKS::tau_f(const BasicField<Scalar, PatchType> &rho_) {
    return muRef / (R * TRef * rho_);
}

/// Initialize
void CDUGKS::initialize() {
    // mesh
    mesh_.output();
    // zone
    forConstRef(zonePatch, mesh_.zones()) {
        auto zonePatchParam = config_.zone(zonePatch.name());
        forConstRef(cellId, zonePatch.group()) {
            auto cell = mesh_.cell(cellId);
            if (cell.rank() != MPI::rank) continue;
            auto rho_ = zonePatchParam.get<Scalar>("density");
            auto U_ = zonePatchParam.get<Vector>("velocity");
            rhoVol.values()[cell.idOnPartition()] = rho_;
            UVol.values()[cell.idOnPartition()] = U_;
        }
    }
    rhoVolOld = rhoVol;
    UVolOld = UVol;
    output();
    // distribution function
    forAll(gVol_, dvId) {
        auto &g = gVol_[dvId];
        const auto &dv = DV_.cell(dvId);
        g = gMaxwell(rhoVol, UVol, dv);
    }
    logger.info << name_ << " - Time: " << time_.name() << std::endl;
    logger.info << "\tInitialize finished." << std::endl;
}

/// Solution
bool CDUGKS::solution() {
    time_.runStep();
    logger.info << name_ << " - Time: " << time_.name() << std::endl;
    dt_ = time_.deltaT();

    updateGbpSurf();

    updateGSurf();

    updateBC();

    updateFVM();

    time_.update();
    return not time_.isStoppable();
}

void CDUGKS::updateGbpSurf() {
    forAll(DV_.cells(), dvId) {
        const auto &dv = DV_.cell(dvId);
        auto &gVol = gVol_[dvId];
        auto &gSurf = gSurf_[dvId];
        auto gVolAdj = fvm::processorCommAdjData(gVol);
        auto gradGVol = fvm::grad(gVol, gVolAdj, fvm::LEAST_SQUARE);
        forConstRef(faceId, gSurf.index()) {
            const auto &face = mesh_.face(faceId);
            const auto nv = face.Snv() / face.S();
            if (face.rank() != MPI::rank) continue;
            const auto &owner = mesh_.cell(face.owner());
            const auto &neighbor = mesh_.cell(face.neighbor());
            if ((nv * dv.C()) > VSMALL) {
                // from owner
                gSurf.values()[face.idOnPartition()] = gVolAdj[face.owner()]
                                                       + gradGVol.values()[owner.idOnPartition()]
                                                         * (face.C() - owner.C() - 0.5 * dt_ * dv.C());
            } else if (nv * dv.C() < -VSMALL) {
                // from neighbor
                gSurf.values()[face.idOnPartition()] = gVolAdj[face.neighbor()]
                                                       + gradGVol.values()[neighbor.idOnPartition()]
                                                         * (face.C() - neighbor.C() - 0.5 * dt_ * dv.C());
            } else {
                // between
                gSurf.values()[face.idOnPartition()] = 0.5 * (gVolAdj[face.owner()] + gVolAdj[face.neighbor()]
                                                              + (gradGVol.values()[owner.idOnPartition()]
                                                                 * (face.C() - owner.C() - 0.5 * dt_ * dv.C())
                                                                 + gradGVol.values()[neighbor.idOnPartition()]
                                                                   * (face.C() - neighbor.C() - 0.5 * dt_ * dv.C())));
            }
        }
    }
}

void CDUGKS::updateGSurf() {
    // Macro var
    rhoSurf = 0.0;
    surfVectorField rhoUSurf(rhoSurf.mesh(), rhoSurf.index());
    forAll(DV_.cells(), dvId) {
        const auto &dv = DV_.cell(dvId);
        const auto &gSurf = gSurf_[dvId];
        rhoSurf += dv.V() * gSurf;
        rhoUSurf += dv.V() * dv.C() * gSurf;
    }
    USurf = rhoUSurf / rhoSurf;
    // recover gSurf
    auto tau = tau_f(rhoSurf);
    auto half_dt = 0.5 * dt_;
    auto factor = half_dt / (2.0 * tau + half_dt);
    forAll(DV_.cells(), dvId) {
        const auto &dv = DV_.cell(dvId);
        auto &gSurf = gSurf_[dvId];
        const auto gEq = gMaxwell(rhoSurf, USurf, dv);
        gSurf = factor * gEq + (1.0 - factor) * gSurf;
    }
}

void CDUGKS::updateBC() {
    forConstRef(markPatch, mesh_.marks()) {
        const auto &markParm = config_.mark(markPatch.name());
        logger.debug << "BC - name: " << markPatch.name() << std::endl;
        switch (Boundary::TypeName.at(markParm.type())) {
            case Boundary::Wall: {
                logger.debug << "BC - wall" << std::endl;
                forConstRef(faceId, markPatch.group()) {
                    const auto &face = mesh_.face(faceId);
                    if (face.rank() != MPI::rank) continue;
                    const auto UW = markParm.get<Vector>("velocity");
                    auto rhoW = 0.0, rhoW0 = 0.0;
                    forAll(DV_.cells(), dvId) {
                        auto &gSurf = gSurf_[dvId].values()[face.idOnPartition()];
                        const auto &dv = DV_.cell(dvId);
                        const auto nv = face.Snv() / face.S();
                        const auto nvk = nv * dv.C();
                        if (nvk > VSMALL) {
                            // from owner
                            rhoW += nvk * gSurf;
                        } else {
                            rhoW0 -= nvk * gMaxwell(1.0, UW, dv);
                        }
                    }
                    rhoW = rhoW / rhoW0;
                    logger.debug << "BC - wall rho=" << rhoW << std::endl;
                    forAll(DV_.cells(), dvId) {
                        auto &gSurf = gSurf_[dvId];
                        const auto &dv = DV_.cell(dvId);
                        const auto nv = face.Snv() / face.S();
                        const auto nvk = nv * dv.C();
                        if (nvk <= VSMALL) {
                            gSurf.values()[face.idOnPartition()] = gMaxwell(rhoW, UW, dv);
                        }
                    }
                }
                break;
            }
            default:
            case Boundary::Processor:
            case Boundary::FluidInterior:
                break;
        }
    }
}

void CDUGKS::updateFVM() {
    logger.debug << "update FVM" << std::endl;
    // Flux
    volScalarField fluxM0(rhoVol.mesh(), rhoVol.index());
    volVectorField fluxM1(rhoVol.mesh(), rhoVol.index());
    forConstRef(cellId, rhoVol.index()) {
        const auto &cell = mesh_.cell(cellId);
        fluxM0.values()[cell.idOnPartition()] = 0.0;
        fluxM1.values()[cell.idOnPartition()] = {0.0, 0.0, 0.0};
        List<List<Scalar>> gSurfAdj(DV_.cellNum(), List<Scalar>(mesh_.faceNum()));
        forAll(DV_.cells(), dvId) {
            gSurfAdj[dvId] = fvm::processorCommAdjData(gSurf_[dvId]);
        }
        forAll(DV_.cells(), dvId) {
            auto &fluxGVol = fluxGVol_[dvId];
            const auto &dv = DV_.cell(dvId);
            auto fluxG = 0.0;
            forConstRef(faceId, cell.faces()) {
                const auto &face = mesh_.face(faceId);
                fluxG += (face.Snv() * dv.C()) * gSurfAdj[dvId][face.id()];
            }
            fluxGVol.values()[cell.idOnPartition()] = fluxG;
            fluxM0.values()[cell.idOnPartition()] += dv.V() * fluxG;
            fluxM1.values()[cell.idOnPartition()] += dv.V() * fluxG * dv.C();
        }
        std::cout << "rank-" << MPI::rank << "flux - cell " << cell.id() << std::endl;
    }
    logger.debug << "flux" << std::endl;
    rhoVolOld = rhoVol;
    UVolOld = UVol;
    volScalarField dt_v(rhoVol.mesh(), rhoVol.index());
    auto rhoUVolOld = rhoVolOld * UVolOld;
    // Marco var
    forConstRef(cellId, rhoVol.index()) {
        const auto &cell = mesh_.cell(cellId);
        dt_v.values()[cell.idOnPartition()] = dt_ / cell.V();
    }
    rhoVol = rhoVolOld - dt_v * fluxM0;
    auto rhoUVol = rhoUVolOld - dt_v * fluxM1;
    UVol = rhoUVol / rhoVol;
    logger.debug << "macro" << std::endl;
    // Evolution
    auto tauOld = tau_f(rhoVolOld);
    auto tau = tau_f(rhoVol);
    auto factor1 = dt_ / (dt_ - 4.0 * tauOld);
    auto factor2 = tau / (tau + 0.5 * dt_);
    auto factor3 = dt_ / (4.0 * tau);
    forAll(DV_.cells(), dvId) {
        const auto &dv = DV_.cell(dvId);
        auto &g = gVol_[dvId];
        auto &fluxG = fluxGVol_[dvId];
        auto gEqOld = gMaxwell(rhoVolOld, UVolOld, dv);
        auto gEq = gMaxwell(rhoVol, UVol, dv);
        auto gOld = factor1 * gEqOld + (1.0 - factor1) * gVol_[dvId];

        g = factor2 * (gOld + dt_ * 0.5 * (gEq / tau + (gEqOld - gOld) / tauOld) - dt_v * fluxG);

        g = factor3 * gEq + (1.0 - factor3) * g;    // g -> gBarPlus
    }
    logger.debug << "evolution" << std::endl;
}

/// Output
void CDUGKS::output() {
    rhoVol.output("rho");
    UVol.output("U");
}
