#include "solver.h"


using Scheme = GKS;
using SCell = Scheme::Cell;
using SFace = Scheme::Face;


/// Solver Constructor
Scheme::GKS(ConfigReader &_config, ArgParser &_parser) : BasicSolver(_config, _parser) {
    /// Read config
    is_prandtle_fix = config.get<bool>("is_prandtle_fix");
    Ma = config.get<double>("Ma");
    gks_solver_key = gks_solver_map.at(config.get<std::string>("gks_solver"));

    mesh.load("./mesh/" + config.phy_mesh);
    mesh.build();
    config.set_mesh_mark(mesh);

    D = mesh.dimension();
}

void Scheme::info() const {
    mesh.info();
}

/// Physical Formula
void Scheme::Prim_to_Cons(PhyVar &_var) const {
    _var.momentum = _var.density * _var.velocity;
    _var.energy = _var.density * _var.velocity * _var.velocity / 2.0 + _var.pressure / (gamma - 1.0);
}

void Scheme::Cons_to_Prim(PhyVar &_var) const {
    _var.velocity = _var.momentum / _var.density;
    _var.pressure = (_var.energy - _var.momentum * _var.momentum / _var.density / 2.0) * (gamma - 1.0);
}

void Scheme::Cons_to_Char(PhyVar &_var) const {
    /// for 2D case
    double over_c = 1.0 / sqrt(gamma * _var.pressure / _var.density);

}

/// Scheme Cell
SCell::Cell(MESH::Cell<int> &cell, Scheme &_solver) : mesh_cell_ptr(&cell), solver(_solver) {
    /// least square
    lsp = generate_least_square(mesh_cell_ptr->key, solver.mesh);
}

void SCell::update() {

}

/// Scheme Face
SFace::Face(MESH::Face<int> &face, Scheme &_solver) : mesh_face_ptr(&face), solver(_solver) {
    is_reduce_order = false;
}

void SFace::reconstruct() {
    /// compute gradient
}

void SFace::calculate_flux() {
    solver.Cons_to_Prim(pv_left);
    solver.Cons_to_Prim(pv_right);
    solver.Cons_to_Prim(pv_center);
    /// correct lambda
    double lambda_left = 0.0;
    double lambda_right = 0.0;
    double lambda_center = 0.0;

    double convar_left[5] = {
            pv_left.density, pv_left.momentum.x, pv_left.momentum.y, pv_right.momentum.z, pv_left.energy
    };
    double convar_right[5] = {
            pv_right.density, pv_right.momentum.x, pv_right.momentum.y, pv_right.momentum.z, pv_right.energy
    };
    double convar0[5] = {
            pv_center.density, pv_center.momentum.x, pv_center.momentum.y, pv_center.momentum.z, pv_center.energy
    };
    /// MMDF

    ///
    if (solver.gks_solver_key == solver.kfvs1st) {
        ;
    }
}
