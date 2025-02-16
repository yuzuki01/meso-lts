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
          rho(mesh_),
          U(mesh_),
          gVol_(DV_.cellNum(), volScalarField(mesh_)),
          gSurf_(DV_.cellNum(), surfScalarField(mesh_)) {
    mesh_.info();
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
            rho.values()[cell.idOnPartition()] = rho_;
            U.values()[cell.idOnPartition()] = U_;
        }
    }
    forAll(gVol_, dvId) {
        auto &g = gVol_[dvId];
        auto dv = DV_.cell(dvId);
        gEqMaxwell(g, rho, U, dv);
    }
    // reduce
    rho = 0.0;
    volVectorField rhoU(rho.mesh(), {0, 0, 0});
    forAll(gVol_, dvId) {
        auto dv = DV_.cell(dvId);
        logger.debug << "dv cell V=" << dv.V() << "  g=" << gVol_[dvId].values()[0] << std::endl;
        rho += dv.V() * gVol_[dvId];
        rhoU += dv.V() * dv.C() * gVol_[dvId];
    }
    U = rhoU / rho;
    rho.output("rho");
    U.output("U");
}

/// Solution
bool CDUGKS::solution() {
    time_.runStep();
    return time_.isStoppable();
}

/// Output
void CDUGKS::output() {

}

/// CDUGKS
void CDUGKS::gEqMaxwell(MESO::volScalarField &g_, const MESO::volScalarField &rho_, const MESO::volVectorField &U_,
                        const Mesh::Cell &dv_) {
    auto RT2 = 2 * R * TRef;
    auto c = dv_.C() - U_;
    g_ = rho_ / pow((M_PI * RT2), mesh_.dimension() * 0.5) * exp(-(c * c) / RT2);
}
