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
          gVol_(DV_.cellNum(), volScalarField(mesh_)),
          gSurf_(DV_.cellNum(), surfScalarField(mesh_)) {
    mesh_.info();
}

/// CDUGKS
template<Label PatchType>
BasicField<Scalar, PatchType> CDUGKS::gMaxwell(const BasicField<Scalar, PatchType> &rho_,
                                               const BasicField<Vector, PatchType> &U_,
                                               const Mesh::Cell &dv_) {
    auto RT2 = 2 * R * TRef;
    auto c = dv_.C() - U_;
    return rho_ / pow((M_PI * RT2), mesh_.dimension() * 0.5) * exp(-(c * c) / RT2);
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
    // distribution function
    forAll(gVol_, dvId) {
        auto &g = gVol_[dvId];
        auto dv = DV_.cell(dvId);
        g = gMaxwell(rhoVol, UVol, dv);
    }
}

/// Solution
bool CDUGKS::solution() {
    time_.runStep();
    return time_.isStoppable();
}

/// Output
void CDUGKS::output() {

}
