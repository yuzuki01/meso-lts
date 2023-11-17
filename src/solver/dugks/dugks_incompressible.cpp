#include "solver.h"

using Scheme = DUGKS_INCOMPRESSIBLE;
using SCell = Scheme::Cell;
using SFace = Scheme::Face;

/// Solver Constructor
Scheme::DUGKS_INCOMPRESSIBLE(ConfigReader & _config, ArgParser & _parser) : config(_config), parser(_parser) {
    /// Logger
    logger(config.solver);
    /// File
    create_dir("./result/" + config.name);
    /// Parser
    max_step = parser.parse_param<int>("max_step", 1e5);
    save_interval = parser.parse_param<int>("save_interval", 1000);
    step = parser.parse_param<int>("continue_step", 0);
    is_crashed = false;
    /// Read config
    Rho = config.get<double>("density");
    Re = config.get<double>("Re");
    Ma = config.get<double>("Ma");
    R = config.get<double>("R");
    T = config.get<double>("temperature");
    L = config.get<double>("length");
    CFL = config.get<double>("CFL");

    phy_mesh.load("./mesh/" + config.phy_mesh);
    phy_mesh.build();
    config.set_mesh_mark(phy_mesh);
    D = phy_mesh.dimension();
    if (config.dvs_mesh == MESH_GAUSS_HERMIT) {
        int gp;
        gp = config.get<int>("gauss_point");
        dvs_mesh = GENERATOR::gauss_hermit(gp, D, R * T);
    } else {
        dvs_mesh.load("./mesh/" + config.dvs_mesh);
        dvs_mesh.build();
    }
}

void Scheme::info() {
    phy_mesh.info();
    dvs_mesh.info();
}

/// Scheme Cell Constructor
SCell::Cell(MESH::Cell<int> & cell, Scheme & _solver) : mesh_cell_ptr(&cell), solver(_solver) {
    f_t.resize(solver.dvs_mesh.cell_num(), 0.0);
    f_bp.resize(solver.dvs_mesh.cell_num(), 0.0);
    slope_f.resize(solver.dvs_mesh.cell_num(), {0.0, 0.0, 0.0});
    /// least square
    update_least_square();
}

/// Scheme Face Constructor
SFace::Face(MESH::Face<int> & face, Scheme & _solver) : mesh_face_ptr(&face), solver(_solver) {
    f.resize(solver.dvs_mesh.cell_num(), 0.0);
    f_b.resize(solver.dvs_mesh.cell_num(), 0.0);
}

/// Physical Formula
double Scheme::f_maxwell(double density, const Vec3D &particle_velocity, const Vec3D &flow_velocity) const {
    double uu_RT = flow_velocity * flow_velocity / RT;
    double ku_RT = flow_velocity * particle_velocity / RT;
    return density * (1.0 + ku_RT + (ku_RT * ku_RT - uu_RT) / 2.0);
}

double Scheme::f_maxwell(const MacroVars &macro_var, const Vec3D &particle_velocity) const {
    double uu_RT = macro_var.velocity * macro_var.velocity / RT;
    double ku_RT = macro_var.velocity * particle_velocity / RT;
    return macro_var.density * (1.0 + ku_RT + (ku_RT * ku_RT - uu_RT) / 2.0);
}

/// Cell Functions
void SCell::update_least_square() {
    lsp = generate_least_square(mesh_cell_ptr->key, solver.phy_mesh);
}

void SCell::get_f_bp() {
    double C = 3.0 * solver.half_dt / (2.0 * solver.tau + solver.dt);
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        f_bp[p] = (1.0 - C) * f_t[p] + C * solver.f_maxwell(macro_vars, it.position);
    }
}

void SCell::get_grad_f_bp() {
    const int near_num = mesh_cell_ptr->near_cell_key.size();
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        Vec3D Sfr(0.0, 0.0, 0.0);
        for (int j = 0; j < near_num; j++) {
            auto &near_cell = solver.get_cell(mesh_cell_ptr->near_cell_key[j]);
            Sfr += lsp.weight[j] * (near_cell.f_bp[p] - f_bp[p]) * lsp.dr[j];
        }
        slope_f[p] = {lsp.Cx * Sfr, lsp.Cy * Sfr, lsp.Cz * Sfr};
        /// zero gradient
        // slope_f[p] = {0.0, 0.0, 0.0};
    }
}

void SCell::get_macro_var() {
    macro_vars.density = macro_vars.temperature = 0.0;
    macro_vars.velocity = macro_vars.heat_flux = {0.0, 0.0, 0.0};
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        macro_vars.density += it.volume * f_t[p];
        macro_vars.velocity += it.volume * f_t[p] * it.position;
    }
    macro_vars.velocity /= macro_vars.density;
}

void SCell::update_f_t() {
    /// compute flux
    std::vector<double> flux_f(solver.dvs_mesh.cell_num(), 0.0);
    /// list face
    double Adt_V;
    for (auto &face_key : mesh_cell_ptr->face_key) {
        auto &face = solver.get_face(face_key);
        auto &mesh_face = solver.phy_mesh.get_face(face_key);
        auto &nv = (mesh_face.on_cell_key == mesh_cell_ptr->key) ? mesh_face.on_cell_nv : mesh_face.inv_cell_nv;
        Adt_V = mesh_face.area / mesh_cell_ptr->volume * solver.dt;
        for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
            auto &it = solver.dvs_mesh.get_cell(p);
            flux_f[p] += (it.position * nv) * Adt_V * face.f[p];
        }
    }
    /// update f_t
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        f_t[p] = (4.0 * f_bp[p] - f_t[p]) / 3.0 - flux_f[p];
    }
}

/// Face Functions
void SFace::get_f_b() {
    auto &on_cell = solver.get_cell(mesh_face_ptr->on_cell_key);
    auto &inv_cell = solver.get_cell(mesh_face_ptr->inv_cell_key);
    auto &nv = mesh_face_ptr->on_cell_nv;
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        if (it.position * nv >= 0.0) {
            Vec3D dr = mesh_face_ptr->position - on_cell.mesh_cell_ptr->position
                       - solver.half_dt * it.position;
            f_b[p] = on_cell.f_bp[p] + on_cell.slope_f[p] * dr;
        } else {
            Vec3D dr = mesh_face_ptr->position - inv_cell.mesh_cell_ptr->position
                       - solver.half_dt * it.position;
            f_b[p] = inv_cell.f_bp[p] + inv_cell.slope_f[p] * dr;
        }
    }
}

void SFace::get_macro_var() {
    macro_vars.density = macro_vars.temperature = 0.0;
    macro_vars.velocity = macro_vars.heat_flux = {0.0, 0.0, 0.0};
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        macro_vars.density += it.volume * f_b[p];
        macro_vars.velocity += it.volume * f_b[p] * it.position;
    }
    macro_vars.velocity /= macro_vars.density;
}

void SFace::get_f() {
    double C = solver.half_dt / (2.0 * solver.tau + solver.half_dt);
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        f[p] = (1.0 - C) * f_b[p] + C * solver.f_maxwell(macro_vars, it.position);
    }
}

void SFace::do_boundary() {
    switch (mesh_face_ptr->boundary_type) {
        case MESH_BC_ISOTHERMAL_WALL: {
            /// isothermal wall
            double kn, rho_w, rho_w0;
            auto &mark = solver.phy_mesh.get_mark(mesh_face_ptr->mark_key);
            auto &nv = mesh_face_ptr->inv_cell_nv;
            rho_w = rho_w0 = 0.0;
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.get_cell(p);
                kn = nv * it.position;
                if (kn >= 0.0) {
                    rho_w0 += it.volume * kn * solver.f_maxwell(1.0, it.position, mark.velocity);
                } else {
                    rho_w -= it.volume * kn * f[p];
                }
            }
            rho_w /= rho_w0;
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.get_cell(p);
                kn = nv * it.position;
                if (kn >= 0.0) {
                    f[p] = solver.f_maxwell(rho_w, it.position, mark.velocity);
                }
            }
            return;
        }
        case MESH_BC_INTERFACE: {
            /// interface - skip
            return;
        }
        default: {
            solver.logger << "Caught unsupported boundary type: " << mesh_face_ptr->info();
            solver.logger.error();
            throw std::invalid_argument("boundary error");
        }
    }
}

/// Scheme Functions
void Scheme::do_step() {
#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = get_cell(i);
        cell.get_f_bp();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = get_cell(i);
        cell.get_grad_f_bp();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = get_face(i);
        face.get_f_b();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = get_face(i);
        face.get_macro_var();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = get_face(i);
        face.get_f();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = get_face(i);
        face.do_boundary();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = get_cell(i);
        cell.update_f_t();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = get_cell(i);
        cell.get_macro_var();
    }

    /// OpenMP loop end
    step++;
    if (step >= max_step && !is_crashed) {
        continue_to_run = false;
        logger << "reach max_step=" << max_step;
        logger.highlight();
        do_save();
        return;
    }
    if (step % save_interval == 0) do_save();
}

void Scheme::init() {
    /// if no init, then cannot run
    continue_to_run = true;
    /// Compute Params
    RT = R * T;
    tau = (Ma * L) / (Re * sqrt(RT));
    dt = CFL * phy_mesh.min_mesh_size / dvs_mesh.max_discrete_velocity;
    half_dt = dt / 2.0;

    logger << "solver info:";
    logger.note();
    data_double_println({"density", "temperature", "length", "R"},
                        {Rho, T, L, R});
    data_double_println({"tau", "CFL", "min-cell-size", "max-velocity"},
                        {tau, CFL, phy_mesh.min_mesh_size, dvs_mesh.max_discrete_velocity});

    /// Create Scheme Objects
    logger << "Create scheme-objects.";
    logger.note();
    for (auto &it : phy_mesh.CELLS) {
        CELLS.emplace_back(it, *this);
    }
    for (auto &cell : CELLS) {
        cell.update_least_square();
        for (auto &it2 : dvs_mesh.CELLS) {
            cell.f_t[it2.key] = f_maxwell(Rho, it2.position, {0.0, 0.0, 0.0});
            cell.f_bp[it2.key] = 0.0;
            cell.slope_f[it2.key] = {0.0, 0.0, 0.0};
            cell.get_macro_var();
        }
    }
    logger << "    scheme-objects: cell - ok.";
    logger.info();
    for (auto &face : phy_mesh.FACES) {
        FACES.emplace_back(face, *this);
        /// face init
        /*
        auto &face = get_face(face_key);
        for (auto & it2 : dvs_mesh.CELLS) {
            auto &key = it2.first;
            auto &particle = it2.second;
            face.f[key] = 0.0;
            face.f_b[key] = 0.0;
        }
         */
    }
    logger << "    scheme-objects: face - ok.";
    logger.info();

    omp_set_num_threads(config.thread);
    logger << "init - ok.";
    logger.note();
}

void Scheme::do_save() {
    std::stringstream ss;
    ss << "./result/" << config.name << "/step_" << step << ".dat";
    string_vector values;
    if (phy_mesh.dimension() == 2) values = {"density", "velocity-x", "velocity-y"};
    else values = {"density", "velocity-x", "velocity-y", "velocity-z"};
    MeshWriter<int> writer(ss.str(), phy_mesh);
    writer.write_head(values);
    writer.write_node();

    std::vector<double> data(phy_mesh.cell_num());

    /// density
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.macro_vars.density;
    writer.write_data(data);
    data.clear();
    /// velocity-x
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.macro_vars.velocity.x;
    writer.write_data(data);
    data.clear();
    /// velocity-y
    for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.macro_vars.velocity.y;
    writer.write_data(data);
    data.clear();

    /// velocity-x
    if (D == 3) {
        for (auto &it : CELLS) data[it.mesh_cell_ptr->key] = it.macro_vars.velocity.z;
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
