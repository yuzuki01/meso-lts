#include "solver.h"


using Scheme = GKS;
using SCell = Scheme::Cell;
using SFace = Scheme::Face;


/// Solver Constructor
Scheme::GKS(ConfigReader &_config, ArgParser &_parser) : BasicSolver(_config, _parser) {
    /// Read config
    is_prandtle_fix = config.get<bool>("is_prandtle_fix");
    Ma = config.get<double>("Ma");
    Re = config.get<double>("Re");
    Pr = config.get<double>("Pr");
    gks_solver_key = gks_solver_map.at(config.get<std::string>("gks_solver"));
    K = config.get<int>("K");
    stage = config.get<int>("stage");
    CFL = config.get<double>("CFL");
    R = config.get<double>("R");
    T = config.get<double>("temperature");
    Rho = config.get<double>("density");
    L = config.get<double>("length");

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

void SCell::update(int now_stage) {

}

/// Scheme Face
SFace::Face(MESH::Face<int> &face, Scheme &_solver) : mesh_face_ptr(&face), solver(_solver) {
    is_reduce_order = false;
}

void SFace::reconstruct(int now_stage) {
    /// compute gradient

}

/*

void SFace::calculate_flux(int now_stage) {
    solver.Cons_to_Prim(pv_left);
    solver.Cons_to_Prim(pv_right);
    solver.Cons_to_Prim(pv_center);

    double Flux[2][5];//XX
    double convar_left[5] = {
            pv_left.density, pv_left.momentum.x, pv_left.momentum.y, pv_right.momentum.z, pv_left.energy
    };
    double convar_right[5] = {
            pv_right.density, pv_right.momentum.x, pv_right.momentum.y, pv_right.momentum.z, pv_right.energy
    };
    double convar0[5] = {
            pv_center.density, pv_center.momentum.x, pv_center.momentum.y, pv_center.momentum.z, pv_center.energy
    };
    //cout << endl;
    //change conservative variables to rho u lambda
    //
    // @param prim_left <rho_L, u_L, v_L, w_L, lambda_L>
    //
    /// double prim_left[5], prim_right[5], prim0[5];
    PhyVar_lambda prim_left, prim_right, prim0;  ///todo set value
    solver.Convar_to_ULambda(prim_left, pv_left);
    solver.Convar_to_ULambda(prim_right, pv_right);
    /// MMDF
    //prim_left[1], prim_left[2], prim_left[3], means U, V, W, Lambda
    MMDF ml(prim_left[1], prim_left[2], prim_left[3], prim_left[4]);
    MMDF mr(prim_right[1], prim_right[2], prim_right[3], prim_right[4]);

    //// Center_all_collision_multi
    solver.Collision(pv_center, prim_left[0], prim_right[0], ml, mr); // ml, mr, pass the whole class into the function
    double a0[5] = { 1.0, 0.0, 0.0, 0.0, 0.0 };
    double alx_t[5], arx_t[5];
    //w_x
    solver.A_point(alx_t,
                   {grad_left.density.x, grad_left.momentum_x.x, grad_left.momentum_y.x, grad_left.momentum_z.x, grad_left.energy.x},
                   prim_left); // input the slope of macroscopic variables, output the a coefficient
    solver.A_point(arx_t,
                   {grad_right.density.x, grad_right.momentum_x.x, grad_right.momentum_y.x, grad_right.momentum_z.x, grad_right.energy.x},
                   prim_right);
    //
    // @param al0x_t <ax1, ax2, ax3, ax4, ax5>
    // Appendix.A (A.2)
    //
    double al0x_t[5];
    double ar0x_t[5];
    solver.GL_address(0, 0, 0, 0, al0x_t, alx_t, ml); // ml, mr, pass the whole class into the function
    solver.GR_address(0, 0, 0, 0, ar0x_t, arx_t, mr); // ml, mr, pass the whole class into the function
    //GL_address, GR_address, G_address, are calculating the macroscopic variables by the microscopic variables based on moment calculation
    //Note: when the input of moment calculation is a coefficient, the result is the derivative of W; when the input is (0 0 0 0 0), the result is W
    //when the input is (1 0 0 0 0), the result is the relative flux, without considering the integral of time

    //w_y
    double aly_t[5], ary_t[5];
    solver.A_point(aly_t,
                   {grad_left.density.y, grad_left.momentum_x.y, grad_left.momentum_y.y, grad_left.momentum_z.y, grad_left.energy.y},
                   prim_left);
    solver.A_point(ary_t,
                   {grad_right.density.y, grad_right.momentum_x.y, grad_right.momentum_y.y, grad_right.momentum_z.y, grad_right.energy.y},
                   prim_right);
    //The difference of A_point with A is just the input form, matrix or pointer
    //The content, input variables, output variables are the same
    double al0y_t[5];
    double ar0y_t[5];
    solver.GL_address(0, 0, 0, 0, al0y_t, aly_t, ml);
    solver.GR_address(0, 0, 0, 0, ar0y_t, ary_t, mr);

    //w_z
    double alz_t[5], arz_t[5];
    solver.A_point(alz_t,
                   {grad_left.density.z, grad_left.momentum_x.z, grad_left.momentum_y.z, grad_left.momentum_z.z, grad_left.energy.z},
                   prim_left);
    solver.A_point(arz_t,
                   {grad_right.density.z, grad_right.momentum_x.z, grad_right.momentum_y.z, grad_right.momentum_z.z, grad_right.energy.z},
            prim_right);

    double al0z_t[5];
    double ar0z_t[5];
    solver.GL_address(0, 0, 0, 0, al0z_t, alz_t, ml);
    solver.GR_address(0, 0, 0, 0, ar0z_t, arz_t, mr);

    for (int var = 0; var < 5; ++var)
    {
        interface.center.der1x[var] = prim_left[0] * al0x_t[var] + prim_right[0] * ar0x_t[var];

        interface.center.der1y[var] = prim_left[0] * al0y_t[var] + prim_right[0] * ar0y_t[var];
        interface.center.der1z[var] = prim_left[0] * al0z_t[var] + prim_right[0] * ar0z_t[var];
    }

    if (solver.gks_solver_key == solver.kfvs1st) {
        ;
    }
}

/// Scheme Functions
void Scheme::do_step() {
    for (int now_stage = 0; now_stage < stage; now_stage++) {
#pragma omp parallel for
        for (int i = 0; i < int(mesh.face_num()); i++) {
            auto &face = get_face(i);
            face.do_boundary(now_stage);     /// ghost-cell
        }
#pragma omp parallel for
        for (int i = 0; i < int(mesh.face_num()); i++) {
            auto &face = get_face(i);
            face.reconstruct(now_stage);
        }
#pragma omp parallel for
        for (int i = 0; i < int(mesh.face_num()); i++) {
            auto &face = get_face(i);
            face.calculate_flux(now_stage);
        }
#pragma omp parallel for
        for (int i = 0; i < int(mesh.cell_num()); i++) {
            auto &cell = get_cell(i);
            cell.update(now_stage);
        }
    }

    /// OpenMP loop end
    step++;
    simulate_time += dt;

    if (step % 10) {
        logger << "step " << step << " time is " << simulate_time;
        logger.info();
    }
    if (simulate_time - stop_time >= 0.0) {
        continue_to_run = false;
        logger << "reach max_simulate_time=" << simulate_time;
        logger.highlight();
        do_save();
        return;
    }
    if (stop_at_specific_time && simulate_time + dt - stop_time > 0.0) {
        dt = stop_time - simulate_time + 1e-16;
    }
    if (step >= max_step && !is_crashed) {
        continue_to_run = false;
        logger << "reach max_step=" << max_step;
        logger.highlight();
        do_save();
        return;
    }
    if (step % save_interval == 0) do_save();
}

*/

/// Solver functions
void Scheme::init() {
    if (D == 2) gamma = (K + 4.0) / (K + 2.0);
    else gamma = (K + 5.0) / (K + 3.0);
    c1_euler = 0.01;
    c2_euler = 1.0;     // e.g. Hypersonic cylinder c2_euler = 5.0
    
    for (auto &it : mesh.CELLS) {
        CELLS.emplace_back(it, *this);
    }
    for (auto &cell : CELLS) {
        // cell.update_geom();
        cell.pv = {};   //todo initial field
    }
    for (auto &it : mesh.FACES) {
        FACES.emplace_back(it, *this);
    }
    for (auto &face : FACES) {
        //todo generate shadow_cell
        //todo generate gaussian point
    }
}


void Scheme::do_save() {
    std::stringstream ss;
    ss << "./result/" << config.name << "/step_" << step << ".dat";
    string_vector values;
    if (mesh.dimension() == 2) values = {"density", "velocity-x", "velocity-y"};    //todo
    else values = {"density", "velocity-x", "velocity-y", "velocity-z"};
    MeshWriter<MESH::ListMesh> writer(ss.str(), mesh);
    writer.write_head(values);
    writer.write_node();

    std::vector<double> data(mesh.cell_num());

    /// density
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.pv.density;
    writer.write_data(data);
    data.clear();
    /// velocity-x
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.pv.velocity.x;
    writer.write_data(data);
    data.clear();
    /// velocity-y
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.pv.velocity.y;
    writer.write_data(data);
    data.clear();

    /// velocity-x
    if (D == 3) {
        for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.pv.velocity.z;
        writer.write_data(data);
        data.clear();
    }

    writer.write_geom();
    writer.close();

    logger << "Save to file: " << ss.str();
    logger.highlight();
}

SCell &Scheme::get_cell(const int &_key) {
    return CELLS[_key];
}

SFace &Scheme::get_face(const int &_key) {
    return FACES[_key];
}
