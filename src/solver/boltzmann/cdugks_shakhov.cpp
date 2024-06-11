#include "solver/solver.h"


using namespace MESO::Solver;


CDUGKS_SHAKHOV::CDUGKS_SHAKHOV(MESO::ArgParser &parser) : BasicSolver(parser) {
    /// params
    Kn = config.get("Kn", 1.0);
    Pr = config.get("Pr", 0.67);
    Ma = config.get("Ma", 0.1);
    CFL = config.get("CFL", 0.8);
    R = config.get("gas-constant", 0.5);
    K = config.get("gas-k", 0);
    Rho0 = config.get("ref-density", 1.0);
    L0 = config.get("ref-length", 1.0);
    T0 = config.get("ref-temperature", 1.0);
    vhs_omega = config.get("vhs-omega", 0.5);
    vhs_index = config.get("vhs-index", 0.81);
    gradient_switch = config.get("gradient-switch", true);
    if (not gradient_switch) {
        logger.warn << "[Warn] gradient-switch was CLOSED." << std::endl;
    }
    limiter_switch = config.get("limiter-switch", false);
    venkata_k = config.get("venkata-limiter-k", 1.0);

    mesh = MESO::Mesh::load_gambit(config.get<std::string>("mesh-file", "<mesh-file>"));
    mesh.info();
    D = mesh.dimension();
    config.info();
    {
        auto dvs_file = config.get<std::string>("dvs-file", "<mesh-file>");
        if (dvs_file == "Newton-Cotes") {
            auto nc_n = config.get("newton-cotes-n", 4);
            auto nc_m = config.get("newton-cotes-mount", 8);
            auto nc_s = config.get("newton-cotes-scale", 12.0);
            dvs_mesh = MESO::Solver::generate_newton_cotes(mesh.dimension(), nc_n, nc_m, nc_s);
        } else {
            dvs_mesh = MESO::Mesh::load_gambit(dvs_file);
        }
    }
    mpi_task = MPI::DVS_partition(dvs_mesh);
    dvs_mesh.info();

    logger.note << "Solver::CDUGKS loaded." << std::endl;
}

void CDUGKS_SHAKHOV::initial() {
    step = 0;
    solution_time = 0.0;
    gamma = (K + 5.0) / (K + 3.0);
    Cv = R / (gamma - 1.0);
    double RT = R * T0;
    double c_sound = sqrt(gamma * RT);
    miu0 = (Kn * L0) * (15 * Rho0 * sqrt(2.0 * M_PI * RT)) / ((7.0 - 2.0 * vhs_omega) * (5.0 - 2.0 * vhs_omega));
    Re = Ma * c_sound * Rho0 * L0 / miu0;
    dt = CFL * mesh.min_cell_size / dvs_mesh.max_cell_magnitude;
    half_dt = dt * 0.5;

    rho_cell = Field<Scalar>(mesh, cell_field_flag);
    T_cell = Field<Scalar>(mesh, cell_field_flag);
    tau_cell = Field<Scalar>(mesh, cell_field_flag);
    vel_cell = Field<Vector>(mesh, cell_field_flag);
    q_cell = Field<Vector>(mesh, cell_field_flag);

    rho_face = Field<Scalar>(mesh, face_field_flag);
    T_face = Field<Scalar>(mesh, face_field_flag);
    tau_face = Field<Scalar>(mesh, face_field_flag);
    vel_face = Field<Vector>(mesh, face_field_flag);
    q_face = Field<Vector>(mesh, face_field_flag);

    g_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    h_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    g_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));
    h_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));

    /// flux
    flux_g.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    flux_h.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));

    auto m0_local = mesh.zero_scalar_field();
    auto m1_local = mesh.zero_vector_field();
    auto m2_local = mesh.zero_scalar_field();
    auto m3_local = mesh.zero_vector_field();
    for (auto &cell: mesh.cells) {
        auto &group = config.get_cell_group(cell, mesh);
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            auto c = particle.position - group.velocity;
            auto cc = c * c;
            auto kk = particle.position * particle.position;
            auto g = g_maxwell(group.density, group.temperature, cc);
            auto h = h_maxwell(group.density, g);
            g_cell[p][cell.id] = g;
            h_cell[p][cell.id] = h;
            m0_local[cell.id] += particle.volume * g;
            m1_local[cell.id] += particle.volume * g * particle.position;
            m2_local[cell.id] += particle.volume * (kk * g + h);
            m3_local[cell.id] += particle.volume * c * (cc * g + h);
        }
    }
    auto m0 = mesh.zero_scalar_field();
    auto m1 = mesh.zero_vector_field();
    auto m2 = mesh.zero_scalar_field();
    auto m3 = mesh.zero_vector_field();
    MPI::AllReduce(m0_local, m0);
    MPI::AllReduce(m1_local, m1);
    MPI::AllReduce(m2_local, m2);
    MPI::AllReduce(m3_local, m3);
    for (auto &cell : mesh.cells) {
        auto rho = m0[cell.id];
        auto rhoU = m1[cell.id];
        auto rhoE =  0.5 * m2[cell.id];
        auto q = 0.5 * m3[cell.id];
        auto u = rhoU / rho;
        auto T = (rhoE / rho - 0.5 * u * u) / Cv;
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = u;
        T_cell[cell.id] = T;
        tau_cell[cell.id] = tau_f(rho, T);
        q_cell[cell.id] = q;
    }

    rho_cell_n = rho_cell;
    vel_cell_n = vel_cell;
    T_cell_n = T_cell;
    tau_cell_n = tau_cell;
    /// residual
    rho_cell_res = rho_cell;
    vel_cell_res = vel_cell;
    T_cell_res = T_cell;
    q_cell_res = q_cell;

    is_crashed = false;
    logger.note << "Initialization finished." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

inline MESO::Scalar CDUGKS_SHAKHOV::tau_f(double rho, double t) const {
    return miu0 * pow(t / T0, vhs_index) / (rho * R * t);
}

inline MESO::Scalar CDUGKS_SHAKHOV::g_maxwell(double rho, double t, double cc) const {
    if (D == 2) {
        double over_RT2 = 0.5 / (R * t);
        return rho * over_RT2 / M_PI * exp(-cc * over_RT2);
    }
    double RT = R * t;
    return rho / pow(2.0 * M_PI * RT, 1.5) * exp(-cc / (2.0 * RT));
}

inline MESO::Scalar CDUGKS_SHAKHOV::h_maxwell(double t, double gm) const {
    return (K + 3.0 - D) * R * t * gm;
}

inline MESO::Scalar CDUGKS_SHAKHOV::g_shakhov(double rho, double t, double cc, double cq, double gm) const {
    /// g_shakhov = g_maxwell + g_pr
    double RT = R * t;
    double p = rho * RT;
    return (1.0 - Pr) * (cq / (5.0 * p * RT)) * (cc / RT - D - 2.0) * gm + gm;
}

inline MESO::Scalar CDUGKS_SHAKHOV::h_shakhov(double rho, double t, double cc, double cq, double gm) const {
    /// h_shakhov = h_maxwell + h_pr
    double RT = R * t;
    double p = rho * RT;
    double h_m = h_maxwell(t, gm);
    return (1.0 - Pr) * (cq / (5.0 * p)) * ((cc / RT - D) * (K + 3.0 - D) - 2.0 * K) * gm + h_m;
}

void CDUGKS_SHAKHOV::reconstruct() {
    /// get f_bar on face
    for (int p = 0; p < mpi_task.size; ++p) {
        ObjectId dvs_id = p + mpi_task.start;
        auto &particle = dvs_mesh.cells[dvs_id];
        /// cell gradient
        Field <Vector> grad_g = g_cell[p].gradient(gradient_switch);
        Field <Vector> grad_h = h_cell[p].gradient(gradient_switch);
        /// interp to face
        for (auto &face: mesh.faces) {
            Vector &nv = face.normal_vector[0];
            auto &cell = (nv * particle.position >= 0.0) ? mesh.cells[face.cell_id[0]] : mesh.cells[face.cell_id[1]];
            Vector dr_ij = face.position - cell.position;
            double phi_g = 1.0, phi_h = 1.0;
            /// venkata limiter
            if (limiter_switch) {
                phi_g = venkata_limiter(g_cell[p], dr_ij * grad_g[cell.id], face, cell, venkata_k);
                phi_h = venkata_limiter(h_cell[p], dr_ij * grad_h[cell.id], face, cell, venkata_k);
            }
            g_face[p][face.id] = g_cell[p][cell.id]
                                 + (dr_ij - particle.position * half_dt) * (phi_g * grad_g[cell.id]);
            h_face[p][face.id] = h_cell[p][cell.id]
                                 + (dr_ij - particle.position * half_dt) * (phi_h * grad_h[cell.id]);
        }
    }

    {
        /// face macro vars
        for (auto &face: mesh.faces) {
            double m0 = 0.0, m2 = 0.0;
            Vector m1(0.0, 0.0, 0.0);
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto kk = particle.position * particle.position;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m0 += particle.volume * g;
                m1 += particle.volume * g * particle.position;
                m2 += particle.volume * (kk * g + h);
            }
            auto u = m1 / m0;
            auto T = (m2 / m0 - u * u) / (2.0 * Cv);
            rho_face[face.id] = m0;
            vel_face[face.id] = u;
            T_face[face.id] = T;
            tau_face[face.id] = tau_f(m0, T);
        }
        for (auto &face: mesh.faces) {
            Vector m3(0.0, 0.0, 0.0);
            auto u = vel_face[face.id];
            auto tau = tau_face[face.id];
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto c = particle.position - u;
                auto cc = c * c;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m3 += particle.volume * c * (cc * g + h);
            }
            auto cq = tau / (2.0 * tau + half_dt * Pr);
            q_face[face.id] = cq * m3;
        }
    }

    for (auto &face: mesh.faces) {
        /// get original f on face
        auto tau = tau_face[face.id];
        auto c_eq = half_dt / (2.0 * tau + half_dt);
        auto rho = rho_face[face.id];
        auto u = vel_face[face.id];
        auto T = T_face[face.id];
        auto q = q_face[face.id];
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            auto c = particle.position - u;
            auto cc = c * c;
            auto cq = c * q;
            auto g_m = g_maxwell(rho, T, cc);
            auto g_s = g_shakhov(rho, T, cc, cq, g_m);
            auto h_s = h_shakhov(rho, T, cc, cq, g_m);
            g_face[p][face.id] = (1.0 - c_eq) * g_face[p][face.id] + c_eq * g_s;
            h_face[p][face.id] = (1.0 - c_eq) * h_face[p][face.id] + c_eq * h_s;
        }
    }

    {
        for (auto &face: mesh.faces) {
            /// boundary
            auto &mark = config.get_face_group(face, mesh);
            auto &nv = face.normal_vector[1];
            auto &neighbor = mesh.cells[face.cell_id[0]];
            switch (mark.type) {
                case inlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= 0.0) {
                            auto c = particle.position - mark.velocity;
                            auto cc = c * c;
                            double g_m = g_maxwell(mark.density, mark.temperature, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(mark.temperature, g_m);
                        }
                    }
                    break;
                case outlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= 0.0) {
                            auto c = particle.position - vel_cell[neighbor.id];
                            auto cc = c * c;
                            auto rho = rho_cell[neighbor.id];
                            auto T = T_cell[neighbor.id];
                            double g_m = g_maxwell(rho, T, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T, g_m);
                        }
                    }
                    break;
                case wall: {
                    double rho_w, rho_w0;
                    rho_w = rho_w0 = 0.0;
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= 0.0) {
                            auto c = particle.position - mark.velocity;
                            auto cc = c * c;
                            auto T = (mark.temperature <= 0.0) ? T_cell[neighbor.id] : mark.temperature;
                            auto g_m = g_maxwell(1.0, T, cc);
                            rho_w0 += kn * particle.volume * g_m;
                        } else {
                            rho_w -= kn * particle.volume * g_face[p][face.id];
                        }
                    }
                    rho_w /= rho_w0;
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= 0.0) {
                            auto c = particle.position - mark.velocity;
                            auto cc = c * c;
                            auto T = (mark.temperature <= 0.0) ? T_cell[neighbor.id] : mark.temperature;
                            auto g_m = g_maxwell(rho_w, T, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T, g_m);
                        }
                    }
                }
                case fluid_interior:
                default:
                    break;
            }
        }
    }
}

void CDUGKS_SHAKHOV::fvm_update() {
    for (auto &cell: mesh.cells) {
        rho_cell_n[cell.id] = rho_cell[cell.id];
        T_cell_n[cell.id] = T_cell[cell.id];
        tau_cell_n[cell.id] = tau_cell[cell.id];
        vel_cell_n[cell.id] = vel_cell[cell.id];

        double flux_m0 = 0.0, flux_m2 = 0.0;
        Vector flux_m1(0.0, 0.0, 0.0);
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            double flux_g_p = 0.0, flux_h_p = 0.0;
            for (auto face_id: cell.face_id) {
                auto &face = mesh.faces[face_id];
                auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                flux_g_p += (particle.position * nv) * face.area * g_face[p][face.id];
                flux_h_p += (particle.position * nv) * face.area * h_face[p][face.id];
            }
            flux_g[p][cell.id] = flux_g_p;
            flux_h[p][cell.id] = flux_h_p;
            flux_m0 += particle.volume * flux_g_p;
            flux_m1 += particle.volume * flux_g_p * particle.position;
            flux_m2 += particle.volume * ((particle.position * particle.position) * flux_g_p + flux_h_p);
        }
        auto dt_v = dt / cell.volume;
        auto rho_n = rho_cell_n[cell.id];
        auto u_n = vel_cell_n[cell.id];
        auto T_n = T_cell_n[cell.id];
        auto rhoU_n = rho_n * u_n;
        auto rhoE_n = rho_n * (0.5 * (u_n * u_n) + Cv * T_n);
        auto rho = rho_n - dt_v * flux_m0;
        auto rhoU = rhoU_n - dt_v * flux_m1;
        auto u = rhoU / rho;
        auto rhoE = rhoE_n - dt_v * flux_m2 * 0.5;
        auto T = (rhoE / rho - (u * u) * 0.5) / Cv;
        auto tau = tau_f(rho, T);
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = rhoU / rho;
        T_cell[cell.id] = T;
        tau_cell[cell.id] = tau;
    }

    for (auto &cell: mesh.cells) {
        auto dt_v = dt / cell.volume;
        /// tn = n
        auto rho_n = rho_cell_n[cell.id];
        auto u_n = vel_cell_n[cell.id];
        auto T_n = T_cell_n[cell.id];
        auto q_n = q_cell[cell.id];
        auto tau_n = tau_cell_n[cell.id];
        /// tn = n + 1
        auto rho = rho_cell[cell.id];
        auto u = vel_cell[cell.id];
        auto T = T_cell[cell.id];
        auto q = q_cell[cell.id];
        auto tau = tau_cell[cell.id];

        auto C = tau / (tau + half_dt);
        double cm = half_dt / (half_dt - 2.0 * tau_n);
        double c_eq = half_dt / (2.0 * tau);
        Vector m3(0.0, 0.0, 0.0);
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            /// tn = n
            auto c_n = particle.position - u_n;
            auto cc_n = c_n * c_n;
            auto cq_n = c_n * q_n;
            auto g_m_n = g_maxwell(rho_n, T_n, cc_n);
            auto g_s_n = g_shakhov(rho_n, T_n, cc_n, cq_n, g_m_n);
            auto h_s_n = h_shakhov(rho_n, T_n, cc_n, cq_n, g_m_n);
            auto g_n = cm * g_s_n + (1.0 - cm) * g_cell[p][cell.id];
            auto h_n = cm * h_s_n + (1.0 - cm) * h_cell[p][cell.id];
            /// tn = n + 1
            auto c = particle.position - u;
            auto cc = c * c;
            auto cq = c * q;
            auto g_m = g_maxwell(rho, T, cc);
            auto g_s = g_shakhov(rho, T, cc, cq, g_m);
            auto h_s = h_shakhov(rho, T, cc, cq, g_m);
            /// Evolution
            auto g = C * (g_n + half_dt * (g_s / tau + (g_s_n - g_n) / tau_n) - dt_v * flux_g[p][cell.id]);
            auto h = C * (h_n + half_dt * (h_s / tau + (h_s_n - h_n) / tau_n) - dt_v * flux_h[p][cell.id]);
            g_cell[p][cell.id] = (1.0 - c_eq) * g + c_eq * g_s;
            h_cell[p][cell.id] = (1.0 - c_eq) * h + c_eq * h_s;
            /// heat-flux
            m3 += particle.volume * c * (cc * g + h);
        }
        /// update heat-flux at tn=n+1
        q_cell[cell.id] = 0.5 * m3;
    }
}

void CDUGKS_SHAKHOV::do_step() {
    reconstruct();
    fvm_update();
    step++;
    solution_time += dt;
    if (MPI::rank == MPI::main_rank) {
        if (step % residual_interval == 0) {
            auto rho_res = residual(rho_cell_res, rho_cell);
            auto vel_res = residual(vel_cell_res, vel_cell);
            auto T_res = residual(T_cell_res, T_cell);
            auto q_res = residual(q_cell_res, q_cell);
            logger.note << "step: " << step << std::endl;
            if (mesh.dimension() == 2) {
                Utils::print_names_and_values({"Res[rho]", "Res[T]", "Res[u]", "Res[v]", "Res[qx]", "Res[qy]"},
                                              {rho_res, T_res, vel_res.x, vel_res.y, q_res.x, q_res.y});
            } else {
                Utils::print_names_and_values(
                        {"Res[rho]", "Res[T]", "Res[u]", "Res[v]", "Res[w]", "Res[qx]", "Res[qy]", "Res[qz]"},
                        {rho_res, T_res, vel_res.x, vel_res.y, vel_res.z, q_res.x, q_res.y, q_res.z});
            }
        }
    }
}

void CDUGKS_SHAKHOV::output() {
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
        return;
    }
    std::stringstream file_name;
    file_name << "./" << case_name << "/result-" << step;

    auto U = vel_cell.heft(0);
    auto V = vel_cell.heft(1);
    auto W = vel_cell.heft(2);
    auto qx = q_cell.heft(0);
    auto qy = q_cell.heft(1);
    auto qz = q_cell.heft(2);
    if (mesh.dimension() == 2) {
        mesh.output(file_name.str(),
                    {"rho", "T", "U", "V", "qx", "qy"},
                    {&rho_cell, &T_cell, &U, &V, &qx, &qy}, step, solution_time);
    } else {
        mesh.output(file_name.str(),
                    {"rho", "T", "U", "V", "W", "qx", "qy", "qz"},
                    {&rho_cell, &T_cell, &U, &V, &W, &qx, &qy, &qz}, step, solution_time);
    }

    if (output_np) {
        Utils::mkdir(case_name + "/np-data");
        /// output numpy data
        rho_cell.output(case_name + "/np-data/rho");
        T_cell.output(case_name + "/np-data/T");
        vel_cell.output(case_name + "/np-data/vel");
        q_cell.output(case_name + "/np-data/q");
    }
}
