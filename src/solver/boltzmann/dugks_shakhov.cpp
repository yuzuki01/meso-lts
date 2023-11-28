#include "solver.h"

using Scheme = DUGKS_SHAKHOV;
using SCell = Scheme::Cell;
using SFace = Scheme::Face;

/// Solver Constructor
Scheme::DUGKS_SHAKHOV(ConfigReader &_config, ArgParser &_parser) : BasicSolver(_config, _parser) {
    /// Read config
    Kn = config.get<double>("Kn");
    Ma = config.get<double>("Ma");
    Pr = config.get<double>("Pr");
    R = config.get<double>("R");
    Rho0 = config.get<double>("density");
    T0 = config.get<double>("temperature");
    L0 = config.get<double>("length");
    CFL = config.get<double>("CFL");
    vhs_index = config.get<double>("vhs_index");
    K = config.get<int>("K");
    limiter_switch = config.get<bool>("limiter_switch");
    zero_gradient = config.get<bool>("zero_gradient");
    boundary_zero_gradient = config.get<bool>("boundary_zero_gradient");
    if (boundary_zero_gradient) {
        logger << "set boundary-zero-gradient.";
        logger.highlight();
    }
    if (zero_gradient) {
        limiter_switch = false;
        logger << "set zero-gradient.";
        logger.highlight();
    } else if (limiter_switch) {
        ventaka_k = config.get<double>("ventaka_k");
        logger << "limiter switch is open, ventaka_k=" << ventaka_k;
        logger.highlight();
    }
    /// Mesh
    phy_mesh.load("./mesh/" + config.phy_mesh);
    phy_mesh.build();

    D = phy_mesh.dimension();
    dvs_mesh.load("./mesh/" + config.dvs_mesh);
    dvs_mesh.build();
}

void Scheme::info() const {
    phy_mesh.info();
    dvs_mesh.info();
}

/// Scheme Cell Constructor
SCell::Cell(MESH::Cell<int> & cell, Scheme &_solver) : mesh_cell (cell), solver(_solver) {
    const int dvs_num = solver.dvs_mesh.cell_num();
    g_t.resize(dvs_num, 0.0);
    g_bp.resize(dvs_num, 0.0);
    slope_g.resize(dvs_num, {0.0, 0.0, 0.0});
    h_t.resize(dvs_num, 0.0);
    h_bp.resize(dvs_num, 0.0);
    slope_h.resize(dvs_num, {0.0, 0.0, 0.0});
    /// least square
    update_geom();
}

/// Scheme Face Constructor
SFace::Face(MESH::Face<int> & face, Scheme &_solver) : mesh_face (face), solver(_solver) {
    const int dvs_num = solver.dvs_mesh.cell_num();
    g.resize(dvs_num, 0.0);
    g_b.resize(dvs_num, 0.0);
    h.resize(dvs_num, 0.0);
    h_b.resize(dvs_num, 0.0);
}

/// Physical Formula
double Scheme::tau_f(double density, double temperature) const {
    return miu0 * pow(temperature / T0, vhs_index) / (density * R * temperature);
}

double Scheme::g_maxwell(double density, double temperature, const Vec3D &c) const {
    double RT2 = R * temperature * 2.0;
    if (D == 2) {
        return density * exp(-(c * c) / RT2) / (Pi * RT2);
    } else {
        return density * exp(-(c * c) / RT2) / pow(Pi * RT2, D / 2.0);
    }
}

double Scheme::h_maxwell(double g_m, double temperature) const {
    return (K + 3.0 - D) * (R * temperature) * g_m;
}

double Scheme::g_shakhov(double density, double temperature, const Vec3D &c, const Vec3D &heat_flux) const {
    double g_m = g_maxwell(density, temperature, c);
    double cq = c * heat_flux;
    double cc = c * c;
    double RT = R * temperature;
    return g_m + (1.0 - Pr) * (cq / 5.0 / density / RT / RT) * (cc / RT - D - 2.0) * g_m;
}

double Scheme::h_shakhov(double density, double temperature, const Vec3D &c, const Vec3D &heat_flux) const {
    double g_m = g_maxwell(density, temperature, c);
    double cq = c * heat_flux;
    double cc = c * c;
    double RT = R * temperature;
    return h_maxwell(g_m, temperature) +
           (1.0 - Pr) * (cq / 5.0 / density / RT) * ((cc / RT - D) * (K + 3.0 - D) - 2.0 * K) * g_m;
}

/// Cell Functions
void SCell::update_geom() {
    lsp = generate_least_square(mesh_cell.key, solver.phy_mesh);
    if (solver.limiter_switch) {
        double kh = solver.ventaka_k * pow(mesh_cell.volume, solver.D / 3.0);
        ventaka_omega = kh * kh * kh;
    }
}

void SCell::get_f_bp() {
    double C = (solver.half_dt + solver.dt) / (2.0 * solver.tau_f(macro_vars.density, macro_vars.temperature) + solver.dt);
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        Vec3D c = it.position - macro_vars.velocity;
        g_bp[p] = (1.0 - C) * g_t[p] +
                  C * solver.g_shakhov(macro_vars.density, macro_vars.temperature, c, macro_vars.heat_flux);
        h_bp[p] = (1.0 - C) * h_t[p] +
                  C * solver.h_shakhov(macro_vars.density, macro_vars.temperature, c, macro_vars.heat_flux);
    }
}

void SCell::get_grad_f_bp() {
    const int near_num = mesh_cell.near_cell_key.size();
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        if (solver.zero_gradient) {
            slope_g[p] = {0.0, 0.0, 0.0};
            slope_h[p] = {0.0, 0.0, 0.0};
            continue;
        }
        Vec3D Sgr{0.0, 0.0, 0.0}, Shr{0.0, 0.0, 0.0};
        for (int j = 0; j < near_num; j++) {
            auto &near_cell = solver.CELLS[mesh_cell.near_cell_key[j]];
            Sgr += lsp.weight[j] * (near_cell.g_bp[p] - g_bp[p]) * lsp.dr[j];
            Shr += lsp.weight[j] * (near_cell.h_bp[p] - h_bp[p]) * lsp.dr[j];
        }
        slope_g[p] = {lsp.Cx * Sgr, lsp.Cy * Sgr, lsp.Cz * Sgr};
        slope_h[p] = {lsp.Cx * Shr, lsp.Cy * Shr, lsp.Cz * Shr};
        /// zero gradient
        // slope_g[p] = slope_h[p] = {0.0, 0.0, 0.0};
    }
}

void SCell::get_macro_var() {
    macro_vars.density = macro_vars.energy = 0.0;
    macro_vars.velocity = macro_vars.heat_flux = {0.0, 0.0, 0.0};
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        macro_vars.density += it.volume * g_t[p];
        macro_vars.velocity += it.volume * g_t[p] * it.position;
        macro_vars.energy += it.volume * (it.position_square * g_t[p] + h_t[p]);
    }
    macro_vars.velocity /= macro_vars.density;
    macro_vars.energy /= 2.0 * macro_vars.density;
    macro_vars.temperature = (macro_vars.energy - (macro_vars.velocity * macro_vars.velocity) / 2.0) / solver.Cv;
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        Vec3D c = it.position - macro_vars.velocity;
        macro_vars.heat_flux += it.volume * c * ((c * c) * g_t[p] + h_t[p]);
    }
    double tau = solver.tau_f(macro_vars.density, macro_vars.temperature);
    macro_vars.heat_flux = (tau / (2.0 * tau + solver.dt * solver.Pr)) * macro_vars.heat_flux;

    /// check
    if (std::isnan(macro_vars.density) || macro_vars.temperature < 0.0) solver.do_crashed(*this);
}

void SCell::update_f_t() {
    /// compute flux
    const int dvs_num = solver.dvs_mesh.cell_num();
    std::vector<double> flux_g(dvs_num, 0.0), flux_h(dvs_num, 0.0);
    /// list face
    double Adt_V;
    for (auto &face_key : mesh_cell.face_key) {
        auto &face = solver.FACES[face_key];
        auto &nv = (face.mesh_face.on_cell_key == mesh_cell.key) ?
                face.mesh_face.on_cell_nv : face.mesh_face.inv_cell_nv;
        Adt_V = face.mesh_face.area / mesh_cell.volume * solver.dt;
        for (int p = 0; p < dvs_num; p++) {
            auto &it = solver.dvs_mesh.CELLS[p];
            flux_g[p] += (it.position * nv) * Adt_V * face.g[p];
            flux_h[p] += (it.position * nv) * Adt_V * face.h[p];
        }
    }
    /// update f_t
    for (int p = 0; p < dvs_num; p++) {
        g_t[p] = (4.0 * g_bp[p] - g_t[p]) / 3.0 - flux_g[p];
        h_t[p] = (4.0 * h_bp[p] - h_t[p]) / 3.0 - flux_h[p];
    }
}

/// Face Functions
void SFace::get_f_b() {
    auto &nv = mesh_face.on_cell_nv;
    auto &on_cell = solver.CELLS[mesh_face.on_cell_key];
    auto &inv_cell = solver.CELLS[mesh_face.inv_cell_key];
    double phi_g = 1.0, phi_h = 1.0;
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.get_cell(p);
        if (it.position * nv >= 0.0) {
            Vec3D r_ij = mesh_face.position - on_cell.mesh_cell.position;
            double g_i = on_cell.g_bp[p], h_i = on_cell.h_bp[p];
            if (solver.boundary_zero_gradient) {
                g_b[p] = g_i;
                h_b[p] = h_i;
                continue;
            }
            if (solver.limiter_switch) {
                Ventaka limiter_g{g_i, g_i, g_i, r_ij * on_cell.slope_g[p]},
                        limiter_h{h_i, h_i, h_i, r_ij * on_cell.slope_h[p]};
                /// find Wi_max and Wi_min
                for (auto &near_cell_key : on_cell.mesh_cell.near_cell_key) {
                    auto &near_cell = solver.CELLS[near_cell_key];
                    double g_n = near_cell.g_bp[p], h_n = near_cell.h_bp[p];
                    if (g_n > limiter_g.Wi_max) limiter_g.Wi_max = g_n;
                    if (g_n < limiter_g.Wi_min) limiter_g.Wi_min = g_n;
                    if (h_n > limiter_h.Wi_max) limiter_h.Wi_max = h_n;
                    if (h_n < limiter_h.Wi_min) limiter_h.Wi_min = h_n;
                }
                /// compute phi
                phi_g = ventaka_phi(limiter_g, on_cell.ventaka_omega);
                phi_h = ventaka_phi(limiter_h, on_cell.ventaka_omega);
            }
            Vec3D dr = mesh_face.position - on_cell.mesh_cell.position - solver.half_dt * it.position;
            g_b[p] = g_i + phi_g * on_cell.slope_g[p] * dr;
            h_b[p] = h_i + phi_h * on_cell.slope_h[p] * dr;
        } else {
            Vec3D r_ij = mesh_face.position - inv_cell.mesh_cell.position;
            double g_i = inv_cell.g_bp[p], h_i = inv_cell.h_bp[p];
            if (solver.boundary_zero_gradient) {
                g_b[p] = g_i;
                h_b[p] = h_i;
                continue;
            }
            if (solver.limiter_switch) {
                Ventaka limiter_g{g_i, g_i, g_i, r_ij * inv_cell.slope_g[p]},
                        limiter_h{h_i, h_i, h_i, r_ij * inv_cell.slope_h[p]};
                /// find Wi_max and Wi_min
                for (auto &near_cell_key : inv_cell.mesh_cell.near_cell_key) {
                    auto &near_cell = solver.CELLS[near_cell_key];
                    double g_n = near_cell.g_bp[p], h_n = near_cell.h_bp[p];
                    if (g_n > limiter_g.Wi_max) limiter_g.Wi_max = g_n;
                    if (g_n < limiter_g.Wi_min) limiter_g.Wi_min = g_n;
                    if (h_n > limiter_h.Wi_max) limiter_h.Wi_max = h_n;
                    if (h_n < limiter_h.Wi_min) limiter_h.Wi_min = h_n;
                }
                /// compute phi
                phi_g = ventaka_phi(limiter_g, inv_cell.ventaka_omega);
                phi_h = ventaka_phi(limiter_h, inv_cell.ventaka_omega);
            }
            Vec3D dr = mesh_face.position - inv_cell.mesh_cell.position - solver.half_dt * it.position;
            g_b[p] = g_i + phi_g * inv_cell.slope_g[p] * dr;
            h_b[p] = h_i + phi_h * inv_cell.slope_h[p] * dr;
        }
    }
}


void SFace::get_macro_var() {
    macro_vars.density = macro_vars.energy = 0.0;
    macro_vars.velocity = macro_vars.heat_flux = {0.0, 0.0, 0.0};
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        macro_vars.density += it.volume * g_b[p];
        macro_vars.velocity += it.volume * g_b[p] * it.position;
        macro_vars.energy += it.volume * (it.position_square * g_b[p] + h_b[p]);
    }
    macro_vars.velocity /= macro_vars.density;
    macro_vars.energy /= 2.0 * macro_vars.density;
    macro_vars.temperature = (macro_vars.energy - (macro_vars.velocity * macro_vars.velocity) / 2.0) / solver.Cv;
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        Vec3D c = it.position - macro_vars.velocity;
        macro_vars.heat_flux += it.volume * c * ((c * c) * g_b[p] + h_b[p]);
    }
    double tau = solver.tau_f(macro_vars.density, macro_vars.temperature);
    macro_vars.heat_flux = (tau / (2.0 * tau + solver.half_dt * solver.Pr)) * macro_vars.heat_flux;

    /// check
    if (std::isnan(macro_vars.density) || macro_vars.temperature < 0.0) solver.do_crashed(*this);
}

void SFace::get_f() {
    double tau = solver.tau_f(macro_vars.density, macro_vars.temperature);
    double C = solver.half_dt / (2.0 * tau + solver.half_dt);
    for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
        auto &it = solver.dvs_mesh.CELLS[p];
        Vec3D c = it.position - macro_vars.velocity;
        g[p] = (1.0 - C) * g_b[p] +
               C * solver.g_shakhov(macro_vars.density, macro_vars.temperature, c, macro_vars.heat_flux);
        h[p] = (1.0 - C) * h_b[p] +
               C * solver.h_shakhov(macro_vars.density, macro_vars.temperature, c, macro_vars.heat_flux);
    }
}

void SFace::do_boundary() {
    switch (mesh_face.boundary_type) {
        case MESH_BC_INTERFACE: {
            /// interface - skip
            return;
        }
        case MESH_BC_INLET: {
            /// inlet
            double kn;
            auto &mark = solver.phy_mesh.get_mark(mesh_face.mark_key);
            auto &nv = mesh_face.inv_cell_nv;
            auto &on_cell = solver.CELLS[mesh_face.on_cell_key];
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.CELLS[p];
                kn = nv * it.position;
                Vec3D c = it.position - mark.velocity;
                double g_m = solver.g_maxwell(mark.density, mark.temperature, c);
                if (kn >= 0.0) {
                    /// in
                    g[p] = g_m;
                    h[p] = solver.h_maxwell(g_m, mark.temperature);
                }
            }
            return;
        }
        case MESH_BC_OUTLET: {
            /// outlet
            auto &mark = solver.phy_mesh.get_mark(mesh_face.mark_key);
            auto &nv = mesh_face.inv_cell_nv;
            auto &on_cell = solver.CELLS[mesh_face.on_cell_key];
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.CELLS[p];
                if (nv * it.position >= 0.0) {
                    Vec3D c = it.position - on_cell.macro_vars.velocity;
                    double g_m = solver.g_maxwell(on_cell.macro_vars.density, on_cell.macro_vars.temperature, c);
                    g[p] = g_m;
                    h[p] = solver.h_maxwell(g_m, on_cell.macro_vars.temperature);
                }
            }
            return;
        }
        case MESH_BC_ISOTHERMAL_WALL: {
            /// isothermal wall
            double kn, rho_w, rho_w0;
            auto &mark = solver.phy_mesh.get_mark(mesh_face.mark_key);
            auto &nv = mesh_face.inv_cell_nv;
            rho_w = rho_w0 = 0.0;
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.CELLS[p];
                kn = nv * it.position;
                if (kn >= 0.0) {
                    rho_w0 += it.volume * kn * solver.g_maxwell(1.0, mark.temperature, it.position - mark.velocity);
                } else {
                    rho_w -= it.volume * kn * g[p];
                }
            }
            rho_w /= rho_w0;
            for (int p = 0; p < solver.dvs_mesh.cell_num(); p++) {
                auto &it = solver.dvs_mesh.CELLS[p];
                kn = nv * it.position;
                if (kn >= 0.0) {
                    Vec3D c = it.position - mark.velocity;
                    double g_m = solver.g_maxwell(rho_w, mark.temperature, c);
                    g[p] = g_m;
                    h[p] = solver.h_maxwell(g_m, mark.temperature);
                }
            }
            return;
        }
        default: {
            solver.logger << "Caught unsupported boundary type: " << mesh_face.info();
            solver.logger.error();
            throw std::invalid_argument("boundary error");
        }
    }
}


/// Scheme Functions
void Scheme::do_step() {
#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = CELLS[i];
        cell.get_f_bp();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = CELLS[i];
        cell.get_grad_f_bp();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = FACES[i];
        face.get_f_b();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = FACES[i];
        face.get_macro_var();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = FACES[i];
        face.get_f();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.face_num()); i++) {
        auto &face = FACES[i];
        face.do_boundary();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = CELLS[i];
        cell.update_f_t();
    }

#pragma omp parallel for
    for (int i = 0; i < int(phy_mesh.cell_num()); i++) {
        auto &cell = CELLS[i];
        cell.get_macro_var();
    }

    /// OpenMP loop end
    step++;
    if (is_crashed) {
        continue_to_run = false;
        logger << "Boom! step = " << step;
        logger.warn();
        do_save();
        return;
    }
    /// debug
    if (debug_mode && step % 10 == 0) {
        logger << "step = " << step;
        logger.debug();
    }
    if (step >= max_step && !is_crashed) {
        continue_to_run = false;
        logger << "reach max_step=" << max_step;
        logger.highlight();
        do_save();
        return;
    }
    if (step % save_interval == 0) {
        do_residual();
        do_save();
    }
}

void Scheme::init() {
    /// if no init, then cannot run
    continue_to_run = true;
    /// Compute Params
    gamma = (5.0 + K) / (3.0 + K);
    Cv = R / (gamma - 1.0);
    double c_sound = sqrt(gamma * R * T0);
    miu0 = (15.0 * Kn * L0 * Rho0 * sqrt(2.0 * Pi * R * T0)) / (2.0 * (7.0 - 2.0 * vhs_index) * (5.0 - 2.0 * vhs_index));
    Re = Ma * c_sound * Rho0 * L0 / miu0;
    dt = CFL * phy_mesh.min_mesh_size / dvs_mesh.max_discrete_velocity;
    half_dt = dt / 2.0;
    /// Boundary
    if (!config.set_mesh_mark(phy_mesh)) {
        continue_to_run = false;
        return;
    }
    logger << "solver info:";
    logger.note();
    data_double_println({"Re", "Ma", "Kn", "Pr"},
                        {Re, Ma, Kn, Pr});
    data_double_println({"density", "temperature", "length", "R"},
                        {Rho0, T0, L0, R});
    data_double_println({"CFL", "min-cell-size", "max-velocity", "time-step"},
                        {CFL, phy_mesh.min_mesh_size, dvs_mesh.max_discrete_velocity, dt});

    /// Create Scheme Objects
    logger << "Create scheme-objects.";
    logger.note();
    for (auto &it : phy_mesh.CELLS) {
        CELLS.emplace_back(it, *this);
    }
    Vec3D init_velocity{0.0, 0.0, 0.0};
    double init_temperature = T0, init_density = Rho0;
    for (auto &mark : phy_mesh.MARKS) {
        if (MESH::MarkTypeID[mark.type] == MESH_BC_INLET) {
            init_velocity = mark.velocity;
            init_temperature = mark.temperature;
            init_density = mark.density;
            break;
        }
    }
    for (auto &cell : CELLS) {
        for (auto &it2 : dvs_mesh.CELLS) {
            Vec3D c = it2.position - init_velocity;
            double g_m = g_maxwell(init_density, init_temperature, c);
            cell.g_t[it2.key] = g_m;
            cell.slope_g[it2.key] = {0.0, 0.0, 0.0};
            cell.h_t[it2.key] = h_maxwell(g_m, init_temperature);
            cell.slope_h[it2.key] = {0.0, 0.0, 0.0};
        }
        cell.get_macro_var();
        cell.macro_vars.old_density = cell.macro_vars.density;
        cell.macro_vars.old_energy = cell.macro_vars.energy;
        /// least-square & limiter
        cell.update_geom();
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
    if (phy_mesh.dimension() == 2) values = {"density", "temperature", "velocity-x", "velocity-y", "heat_flux-x", "heat_flux-y"};
    else values = {"density", "temperature", "velocity-x", "velocity-y", "velocity-z", "heat_flux-x", "heat_flux-y", "heat_flux-z"};
    MeshWriter<MESH::ListMesh> writer(ss.str(), phy_mesh);
    writer.write_head(values);
    writer.write_node();

    std::vector<double> data(phy_mesh.cell_num());

    /// density
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.density;
    writer.write_data(data);
    data.clear();
    /// temperature
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.temperature;
    writer.write_data(data);
    data.clear();
    /// velocity-x
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.x;
    writer.write_data(data);
    data.clear();
    /// velocity-y
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.y;
    writer.write_data(data);
    data.clear();
    /// velocity-x
    if (D == 3) {
        for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.z;
        writer.write_data(data);
        data.clear();
    }
    /// heat_flux-x
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.heat_flux.x;
    writer.write_data(data);
    data.clear();
    /// heat_flux-y
    for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.heat_flux.y;
    writer.write_data(data);
    data.clear();
    /// heat_flux-x
    if (D == 3) {
        for (auto &it : CELLS) data[it.mesh_cell.key] = it.macro_vars.heat_flux.z;
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

void Scheme::do_crashed(Scheme::Cell &cell) {
    logger << "Caught crashed value: \n"
              "   cell<" << cell.mesh_cell.key << ", pos=" << cell.mesh_cell.position.info() << ">\n"
                                                                                                "   density=" << cell.macro_vars.density << " temperature=" << cell.macro_vars.temperature;
    logger.warn();
    is_crashed = true;
}

void Scheme::do_crashed(Scheme::Face &face) {
    logger << "Caught crashed value (checking before do_boundary()): \n"
              "   face<" << face.mesh_face.key << ", pos=" << face.mesh_face.position.info() << ">\n"
              "      on_cell=" << face.mesh_face.on_cell_key << "  inv_cell=" << face.mesh_face.inv_cell_key << "\n"
              "   density=" << face.macro_vars.density << " temperature=" << face.macro_vars.temperature;
    logger.warn();
    /// is_crashed = true;
}

void Scheme::do_residual() {
    std::vector<double> residual(2, 0.0);
    double sum1, sum2;
    /// density
    sum1 = sum2 = 0.0;
    for (auto &cell : CELLS) {
        double a, b;
        a = cell.macro_vars.density * cell.mesh_cell.volume;
        b = cell.macro_vars.old_density * cell.mesh_cell.volume;
        sum1 += fabs(a * a - b * b);
        sum2 += fabs(b * b);
        cell.macro_vars.old_density = cell.macro_vars.density;
    }
    residual[0] = sqrt(sum1 / sum2);
    /// energy
    sum1 = sum2 = 0.0;
    for (auto &cell : CELLS) {
        double a, b;
        a = cell.macro_vars.energy * cell.mesh_cell.volume;
        b = cell.macro_vars.old_energy * cell.mesh_cell.volume;
        sum1 += fabs(a * a - b * b);
        sum2 += fabs(b * b);
        cell.macro_vars.old_energy = cell.macro_vars.energy;
    }
    residual[1] = sqrt(sum1 / sum2);
    logger << "Residual at step=" << step;
    logger.info();
    data_double_println({"density", "energy"}, residual);
}
