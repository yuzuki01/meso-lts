#include "solver/solver.h"


using namespace MESO::Solver;

/// solver coding
void CDUGKS_SHAKHOV::update_config() {
    config.update_config(false);

    Pr = config.get("Pr", 0.67);
    dt = config.get("timeStep", 1.0e-5);
    R = config.get("gas-constant", 0.5);
    K = config.get("gas-k", 0);
    Rho0 = config.get("ref-density", 1.0);
    L0 = config.get("ref-length", 1.0);
    T0 = config.get("ref-temperature", 1.0);
    miu0 = config.get("ref-viscosity", 1.60e-05);
    vhs_index = config.get("vhs-index", 0.81);
    gradient_switch = config.get("gradient-switch", true);
    limiter_switch = config.get("limiter-switch", false);
    if (not gradient_switch and step == 0) {
        logger.warn << "[Warn] gradient-switch was CLOSED." << std::endl;
        if (not limiter_switch) {
            logger.warn << "[Warn] limiter-switch was CLOSED." << std::endl;
        } else {
            venkata_k = config.get("venkata-limiter-k", 1.0);
        }
    }
    /// calculate params
    gamma = (K + 5.0) / (K + 3.0);
    Cv = (K + 3.0) * R * 0.5;
    half_dt = dt * 0.5;
}

CDUGKS_SHAKHOV::CDUGKS_SHAKHOV(MESO::ArgParser &parser, Config &config) : BasicSolver(parser, config) {
    update_config();
    mesh = MESO::fvmMesh::load_gambit(config.get<String>("mesh-file", "<mesh-file>"),
                                      mesh_scale);
    mesh.info();
    D = mesh.dimension();
    config.info();

    auto dvs_file = config.get<String>("dvs-file", "<mesh-file>");
    auto dvs_type = config.get<String>("dvs-type", "<dvs-type>");
    if (dvs_type == "Newton-Cotes") {
        auto dvs_params = DVS::NewtonCotesParams(mesh.dimension());
        dvs_mesh = DVS::generate_dvs(dvs_file, dvs_params);
    } else if (dvs_type == "half-range-Gauss-Hermite") {
        auto dvs_params = DVS::HalfRangeGHParams(mesh.dimension(), R * T0);
        dvs_mesh = DVS::generate_dvs(dvs_file, dvs_params);
    } else if (dvs_type == "Sparse") {
        dvs_mesh = DVS::read_mesh_file(dvs_file);
    } else {
        dvs_mesh = fvmMesh::load_gambit(dvs_file);
    }
    mpi_task = MPI::DVS_partition(dvs_mesh);
    dvs_mesh.info();

    logger.note << "Solver::CDUGKS loaded." << std::endl;
}

void CDUGKS_SHAKHOV::initial() {
    step = 0;
    solution_time = 0.0;
    gamma = (K + 5.0) / (K + 3.0);
    Cv = (K + 3.0) * R * 0.5;

    rho_cell = mesh.zero_scalar_field(cell_field_flag);
    T_cell = mesh.zero_scalar_field(cell_field_flag);
    vel_cell = mesh.zero_vector_field(cell_field_flag);
    q_cell = mesh.zero_vector_field(cell_field_flag);

    g_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    h_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    g_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));
    h_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));
    /// flux
    flux_g.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    flux_h.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));

    if (config.get("read-np-data", false, false)) {
        /// Init by np-data
        read_np_data<Scalar>(case_name + "/latest/Rho", rho_cell);
        read_np_data<Scalar>(case_name + "/latest/T", T_cell);
        read_np_data<Vector>(case_name + "/latest/U", vel_cell);
        read_np_data<Vector>(case_name + "/latest/q", q_cell);
        for (auto &cell: mesh.cells) {
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto c = particle.position - vel_cell[cell.id];
                auto cc = c * c;
                auto cq = c * q_cell[cell.id];
                auto kk = particle.position * particle.position;
                auto g_m = g_maxwell(rho_cell[cell.id], T_cell[cell.id], cc);
                auto g = g_shakhov(rho_cell[cell.id], T_cell[cell.id], cc, cq, g_m);
                auto h = h_shakhov(rho_cell[cell.id], T_cell[cell.id], cc, cq, g_m);
                g_cell[p][cell.id] = g;
                h_cell[p][cell.id] = h;
            }
        }
    } else {
        auto m0_local = mesh.zero_scalar_field();
        auto m1_local = mesh.zero_vector_field();
        auto m2_local = mesh.zero_scalar_field();
        auto m3_local = mesh.zero_vector_field();
        /// Init by cell group
        for (auto &cell: mesh.cells) {
            auto &group = config.get_cell_group(cell, mesh);

            auto rho_patch_type = group.patch.get_type("density");
            Scalar rho_patch;
            if (rho_patch_type == PatchType::fromFile) {
                rho_patch = group.patch.get_file_scalar("density", cell.id);
            } else {
                rho_patch = group.patch.get_scalar("density");
            }
            auto u_patch_type = group.patch.get_type("velocity");
            Vector u_patch;
            if (u_patch_type == PatchType::fromFile) {
                u_patch = group.patch.get_file_vector("velocity", cell.id);
            } else {
                u_patch = group.patch.get_vector("velocity");
            }
            auto T_patch_type = group.patch.get_type("temperature");
            Scalar T_patch;
            if (T_patch_type == PatchType::fromFile) {
                T_patch = group.patch.get_file_scalar("temperature", cell.id);
            } else {
                T_patch = group.patch.get_scalar("temperature");
            }

            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto c = particle.position - u_patch;
                auto cc = c * c;
                auto kk = particle.position * particle.position;
                auto g = g_maxwell(rho_patch, T_patch, cc);
                auto h = h_maxwell(T_patch, g);
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
        for (auto &cell: mesh.cells) {
            auto rho = m0[cell.id];
            auto rhoU = m1[cell.id];
            auto rhoE = 0.5 * m2[cell.id];
            auto q = 0.5 * m3[cell.id];
            auto u = rhoU / rho;
            auto T = (rhoE / rho - 0.5 * u * u) / Cv;
            rho_cell[cell.id] = rho;
            vel_cell[cell.id] = u;
            T_cell[cell.id] = T;
            q_cell[cell.id] = q;
        }
    }
    /// residual
    rho_cell_res = rho_cell;
    vel_cell_res = vel_cell;
    T_cell_res = T_cell;
    q_cell_res = q_cell;

    run_state = true;
    logger.note << "Initialization finished." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

inline MESO::Scalar CDUGKS_SHAKHOV::tau_f(Scalar rho, Scalar t) const {
    return miu0 * pow(t / T0, vhs_index) / (rho * R * t);
}

inline MESO::Scalar CDUGKS_SHAKHOV::g_maxwell(Scalar rho, Scalar t, Scalar cc) const {
    if (D == 2) {
        Scalar over_RT2 = 0.5 / (R * t);
        return rho * over_RT2 / M_PI * exp(-cc * over_RT2);
    }
    Scalar RT = R * t;
    return rho / pow(2.0 * M_PI * RT, 1.5) * exp(-cc * 0.5 / RT);
}

inline MESO::Scalar CDUGKS_SHAKHOV::h_maxwell(Scalar t, Scalar gm) const {
    return (K + 3.0 - D) * (R * t) * gm;
}

inline MESO::Scalar CDUGKS_SHAKHOV::g_shakhov(Scalar rho, Scalar t, Scalar cc, Scalar cq, Scalar gm) const {
    /// g_shakhov = g_maxwell + g_pr
    Scalar RT = R * t;
    Scalar p = rho * RT;
    return (1.0 - Pr) * (cq / (5.0 * p * RT)) * (cc / RT - D - 2.0) * gm + gm;
}

inline MESO::Scalar CDUGKS_SHAKHOV::h_shakhov(Scalar rho, Scalar t, Scalar cc, Scalar cq, Scalar gm) const {
    /// h_shakhov = h_maxwell + h_pr
    Scalar RT = R * t;
    Scalar p = rho * RT;
    Scalar h_m = h_maxwell(t, gm);
    return (1.0 - Pr) * (cq / (5.0 * p)) * ((cc / RT - D) * (K + 3.0 - D) - 2.0 * K) * gm + h_m;
}

void CDUGKS_SHAKHOV::reconstruct() {
    /// get f_bar on face
    for (int p = 0; p < mpi_task.size; ++p) {
        ObjectId dvs_id = p + mpi_task.start;
        auto &particle = dvs_mesh.cells[dvs_id];
        if (gradient_switch) {
            /// cell gradient
            Field <Vector> grad_g = fvmMesh::grad(g_cell[p]);
            Field <Vector> grad_h = fvmMesh::grad(h_cell[p]);
            /// interp to face
            for (auto &face: mesh.faces) {
                Vector &nv = face.normal_vector[0];
                auto &own = mesh.cells[face.cell_id[0]];
                auto &nei = mesh.cells[face.cell_id[1]];
                auto nv_xi = nv * particle.position;
                if (nv_xi >= VSMALL) {
                    // from own
                    Vector dr_ij = face.position - own.position;
                    /// venkata-limiter
                    Scalar phi_g = 1.0, phi_h = 1.0;
                    if (limiter_switch) {
                        phi_g = venkata_limiter(g_cell[p], dr_ij * grad_g[own.id], own, venkata_k, 1.0e-10);
                        phi_h = venkata_limiter(h_cell[p], dr_ij * grad_h[own.id], own, venkata_k, 1.0e-10);
                    }
                    g_face[p][face.id] = g_cell[p][own.id]
                                         + phi_g * ((dr_ij - particle.position * half_dt) * grad_g[own.id]);
                    h_face[p][face.id] = h_cell[p][own.id]
                                         + phi_h * ((dr_ij - particle.position * half_dt) * grad_h[own.id]);
                } else if (nv_xi < -VSMALL) {
                    // from nei
                    Vector dr_ij = face.position - nei.position;
                    /// venkata-limiter
                    Scalar phi_g = 1.0, phi_h = 1.0;
                    if (limiter_switch) {
                        phi_g = venkata_limiter(g_cell[p], dr_ij * grad_g[nei.id], nei, venkata_k, VSMALL);
                        phi_h = venkata_limiter(h_cell[p], dr_ij * grad_h[nei.id], nei, venkata_k, VSMALL);
                    }
                    g_face[p][face.id] = g_cell[p][nei.id]
                                         + phi_g * ((dr_ij - particle.position * half_dt) * grad_g[nei.id]);
                    h_face[p][face.id] = h_cell[p][nei.id]
                                         + phi_h * ((dr_ij - particle.position * half_dt) * grad_h[nei.id]);
                } else {
                    // between
                    Vector dr_own = face.position - own.position;
                    Vector dr_nei = face.position - nei.position;
                    /// venkata-limiter
                    Scalar phi_g_own = 1.0, phi_h_own = 1.0;
                    Scalar phi_g_nei = 1.0, phi_h_nei = 1.0;
                    if (limiter_switch) {
                        phi_g_own = venkata_limiter(g_cell[p], dr_own * grad_g[own.id], own, venkata_k, VSMALL);
                        phi_h_own = venkata_limiter(h_cell[p], dr_own * grad_h[own.id], own, venkata_k, VSMALL);
                        phi_g_nei = venkata_limiter(g_cell[p], dr_nei * grad_g[nei.id], nei, venkata_k, VSMALL);
                        phi_h_nei = venkata_limiter(h_cell[p], dr_nei * grad_h[nei.id], nei, venkata_k, VSMALL);
                    }
                    g_face[p][face.id] = 0.5 * (g_cell[p][own.id] + g_cell[p][nei.id]
                                                + phi_g_own * ((dr_own - particle.position * half_dt) * grad_g[own.id])
                                                + phi_g_nei * ((dr_nei - particle.position * half_dt) * grad_g[nei.id])
                    );
                    h_face[p][face.id] = 0.5 * (h_cell[p][own.id] + h_cell[p][nei.id]
                                                + phi_h_own * ((dr_own - particle.position * half_dt) * grad_h[own.id])
                                                + phi_h_nei * ((dr_nei - particle.position * half_dt) * grad_h[nei.id])
                    );
                }
            }
        } else {
            /// zeroGradient to face
            for (auto &face: mesh.faces) {
                auto &own = mesh.cells[face.cell_id[0]];
                auto &nei = mesh.cells[face.cell_id[1]];
                Vector &nv = face.normal_vector[0];
                auto nv_xi = nv * particle.position;
                if (nv_xi >= VSMALL) {
                    // from own
                    g_face[p][face.id] = g_cell[p][own.id];
                    h_face[p][face.id] = h_cell[p][own.id];
                } else if (nv_xi < -VSMALL) {
                    // from nei
                    g_face[p][face.id] = g_cell[p][nei.id];
                    h_face[p][face.id] = h_cell[p][nei.id];
                } else {
                    // between
                    g_face[p][face.id] = 0.5 * (g_cell[p][own.id] + g_cell[p][nei.id]);
                    h_face[p][face.id] = 0.5 * (h_cell[p][own.id] + h_cell[p][nei.id]);
                }
            }
        }
    }
    {
        auto m0_local = mesh.zero_scalar_field(face_field_flag);
        auto m1_local = mesh.zero_vector_field(face_field_flag);
        auto m2_local = mesh.zero_scalar_field(face_field_flag);
        /// face macro vars
        for (auto &face: mesh.faces) {
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto kk = particle.position * particle.position;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m0_local[face.id] += particle.volume * g;
                m1_local[face.id] += particle.volume * g * particle.position;
                m2_local[face.id] += particle.volume * (kk * g + h);
            }
        }
        auto m0 = mesh.zero_scalar_field(face_field_flag);
        auto m1 = mesh.zero_vector_field(face_field_flag);
        auto m2 = mesh.zero_scalar_field(face_field_flag);
        MPI::AllReduce(m0_local, m0);
        MPI::AllReduce(m1_local, m1);
        MPI::AllReduce(m2_local, m2);
        auto m3_local = mesh.zero_vector_field(face_field_flag);
        for (auto &face: mesh.faces) {
            auto rho = m0[face.id];
            auto rhoU = m1[face.id];
            auto u = rhoU / rho;
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto c = particle.position - u;
                auto cc = c * c;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m3_local[face.id] += particle.volume * c * (cc * g + h);
            }
        }
        auto m3 = mesh.zero_vector_field(face_field_flag);
        MPI::AllReduce(m3_local, m3);
        /// get original f on face
        for (auto &face: mesh.faces) {
            auto rho = m0[face.id];
            auto rhoU = m1[face.id];
            auto rhoE = 0.5 * m2[face.id];
            auto u = rhoU / rho;
            auto T = (rhoE / rho - 0.5 * u * u) / Cv;
            auto tau = tau_f(rho, T);
            auto q = (tau / (2.0 * tau + half_dt * Pr)) * m3[face.id];
            auto C_s = half_dt / (2.0 * tau + half_dt);
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto c = particle.position - u;
                auto cc = c * c;
                auto cq = c * q;
                auto g_m = g_maxwell(rho, T, cc);
                auto g_s = g_shakhov(rho, T, cc, cq, g_m);
                auto h_s = h_shakhov(rho, T, cc, cq, g_m);
                g_face[p][face.id] = (1.0 - C_s) * g_face[p][face.id] + C_s * g_s;
                h_face[p][face.id] = (1.0 - C_s) * h_face[p][face.id] + C_s * h_s;
            }
        }
    }
    {
        Map <ObjectId> wall_id_map;
        List <Scalar> wall_rho_local, wall_rho0_local;

        /// boundary
        for (auto &face: mesh.faces) {
            auto &mark = config.get_face_group(face, mesh);
            auto &nv = face.normal_vector[1];
            auto &neighbor = mesh.cells[face.cell_id[0]];
            switch (mark.type) {
                case BoundaryType::pressure_inlet: {
                    auto T_m = T_cell[neighbor.id];
                    auto u_m = vel_cell[neighbor.id];
                    auto u_in = (u_m * nv) * nv;
                    auto T_in = mark.patch.get_scalar("temperature");
                    auto p_in = Boundary::Compressible::solve_pressure(
                            mark.patch.get_scalar("pressure"), R, T_m, u_in, gamma);
                    auto rho_in = p_in / (R * T_in);
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            auto c = particle.position - u_in;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_in, T_in, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T_in, g_m);
                        }
                    }
                }
                    break;
                case BoundaryType::farfield_inlet: {
                    auto rho_patch = mark.patch.get_scalar("density");
                    auto u_patch = mark.patch.get_vector("velocity");
                    auto T_patch = mark.patch.get_scalar("temperature");
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_patch, T_patch, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T_patch, g_m);
                        }
                    }
                }
                    break;
                case BoundaryType::freestream_inlet: {
                    auto p_patch = mark.patch.get_scalar("pressure");
                    auto u_patch = mark.patch.get_vector("velocity");
                    auto T_patch = mark.patch.get_scalar("temperature");
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto rho_patch = p / (R * T_patch);
                        if (particle.position * nv >= VSMALL) {
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_patch, T_patch, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T_patch, g_m);
                        }
                    }
                }
                    break;
                case BoundaryType::pressure_outlet: {
                    auto p_patch = mark.patch.get_scalar("pressure");
                    auto u_m = vel_cell[neighbor.id];
                    auto T_m = T_cell[neighbor.id];
                    auto p_e = Boundary::Compressible::solve_pressure(p_patch, R,
                                                                      T_m, u_m, gamma);
                    auto rho_e = p_e / (R * T_m);
                    auto nv_e = -nv;
                    auto u_e = (u_m * nv_e) * nv_e;
                    auto T_e = mark.patch.get_scalar("temperature");
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            auto c = particle.position - u_e;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_e, T_e, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T_e, g_m);
                        }
                    }
                }
                    break;
                case BoundaryType::farfield_outlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            auto c = particle.position - vel_cell[neighbor.id];
                            auto cc = c * c;
                            auto rho = rho_cell[neighbor.id];
                            auto T = T_cell[neighbor.id];
                            auto g_m = g_maxwell(rho, T, cc);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_maxwell(T, g_m);
                        }
                    }
                    break;
                case BoundaryType::slip_wall:
                case BoundaryType::wall: {
                    Scalar rho_w, rho_w0;
                    rho_w = rho_w0 = 0.0;
                    Vector u_w = mark.patch.get_vector("velocity");
                    auto T_patch_type = mark.patch.get_type("temperature");
                    Scalar T_w = 0.0;
                    if (T_patch_type == PatchType::zeroGradient) {
                        T_w = T_cell[neighbor.id];
                    } else if (T_patch_type == PatchType::fixedValue) {
                        T_w = mark.patch.get_scalar("temperature");
                    }
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            auto c = particle.position - u_w;
                            auto cc = c * c;
                            auto g_m = g_maxwell(1.0, T_w, cc);
                            rho_w0 += kn * particle.volume * g_m;
                        } else {
                            rho_w -= kn * particle.volume * g_face[p][face.id];
                        }
                    }
                    wall_id_map[face.id] = int(wall_rho_local.size());
                    wall_rho_local.push_back(rho_w);
                    wall_rho0_local.push_back(rho_w0);
                }
                    break;
                case BoundaryType::symmetry: {
                    auto direction = mark.patch.get_int("direction");
                    List <Scalar> g_all(dvs_mesh.NCELL), h_all(dvs_mesh.NCELL);
                    MPI::GatherFieldList(g_face, g_all, face.id);
                    MPI::GatherFieldList(h_face, h_all, face.id);
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            g_face[p][face.id] = g_all[particle.symmetry.id[direction]];
                            h_face[p][face.id] = h_all[particle.symmetry.id[direction]];
                        } else {
                            if (kn < -VSMALL) continue;
                            g_face[p][face.id] = 0.5 * (g_face[p][face.id] + g_all[particle.symmetry.id[direction]]);
                            h_face[p][face.id] = 0.5 * (h_face[p][face.id] + h_all[particle.symmetry.id[direction]]);
                        }
                    }
                }
                    break;
                case BoundaryType::fluid_interior:
                default:
                    break;
            }
        }
        List <Scalar> wall_rho_global;
        List <Scalar> wall_rho0_global;
        MPI::AllReduce(wall_rho_local, wall_rho_global);
        MPI::AllReduce(wall_rho0_local, wall_rho0_global);

        for (auto &face: mesh.faces) {
            auto &mark = config.get_face_group(face, mesh);
            switch (mark.type) {
                case BoundaryType::wall: {
                    auto u_patch = mark.patch.get_vector("velocity");
                    auto T_patch = mark.patch.get_scalar("temperature");

                    auto &nv = face.normal_vector[1];
                    int wall_list_id = wall_id_map[face.id];
                    auto rho_w = wall_rho_global[wall_list_id] / wall_rho0_global[wall_list_id];
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_w, T_patch, cc);
                            auto h_m = h_maxwell(T_patch, g_m);
                            g_face[p][face.id] = g_m;
                            h_face[p][face.id] = h_m;
                        } else {
                            if (kn < -VSMALL) continue;
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_w, T_patch, cc);
                            auto h_m = h_maxwell(T_patch, g_m);
                            g_face[p][face.id] = 0.5 * (g_m + g_face[p][face.id]);
                            h_face[p][face.id] = 0.5 * (h_m + h_face[p][face.id]);
                        }
                    }
                }
                    break;
                case BoundaryType::slip_wall: {
                    auto u_patch = mark.patch.get_vector("velocity");
                    auto T_patch = mark.patch.get_scalar("temperature");
                    auto direction = mark.patch.get_int("direction");
                    auto alpha_u = mark.patch.get_scalar("slip-wall-alpha-u");
                    auto alpha_T = mark.patch.get_scalar("slip-wall-alpha-T");

                    List <Scalar> g_all(dvs_mesh.NCELL), h_all(dvs_mesh.NCELL);
                    MPI::GatherFieldList(g_face, g_all, face.id);
                    MPI::GatherFieldList(h_face, h_all, face.id);
                    auto &nv = face.normal_vector[1];
                    int wall_list_id = wall_id_map[face.id];
                    auto rho_w = wall_rho_global[wall_list_id] / wall_rho0_global[wall_list_id];
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        auto kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_w, T_patch, cc);
                            auto h_m = h_maxwell(T_patch, g_m);
                            auto g_r = g_all[particle.symmetry.id[direction]];
                            auto h_r = h_all[particle.symmetry.id[direction]];
                            g_face[p][face.id] = (1.0 - alpha_u) * g_m + alpha_u * g_r;
                            h_face[p][face.id] = (1.0 - alpha_T) * h_m + alpha_T * h_r;
                        } else {
                            if (kn < -VSMALL) continue;
                            auto c = particle.position - u_patch;
                            auto cc = c * c;
                            auto g_m = g_maxwell(rho_w, T_patch, cc);
                            auto h_m = h_maxwell(T_patch, g_m);
                            auto g_r = g_all[particle.symmetry.id[direction]];
                            auto h_r = h_all[particle.symmetry.id[direction]];
                            g_face[p][face.id] = 0.5 * ((1.0 - alpha_u) * g_m + alpha_u * g_r
                                                        + g_face[p][face.id]);
                            h_face[p][face.id] = 0.5 * ((1.0 - alpha_T) * h_m + alpha_T * h_r
                                                        + h_face[p][face.id]);
                        }
                    }
                }
                    break;
                default:
                    break;
            }
        }
    }
}

void CDUGKS_SHAKHOV::fvm_update() {
    /// Flux
    auto flux_m0_local = mesh.zero_scalar_field();
    auto flux_m1_local = mesh.zero_vector_field();
    auto flux_m2_local = mesh.zero_scalar_field();

    for (auto &cell: mesh.cells) {
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            Scalar flux_g_tmp = 0.0, flux_h_tmp = 0.0;
            for (auto face_id: cell.face_id) {
                auto &face = mesh.faces[face_id];
                auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                auto knA = (particle.position * nv) * face.area;
                flux_g_tmp += knA * g_face[p][face.id];
                flux_h_tmp += knA * h_face[p][face.id];
            }
            flux_g[p][cell.id] = flux_g_tmp;
            flux_h[p][cell.id] = flux_h_tmp;
            flux_m0_local[cell.id] += particle.volume * flux_g_tmp;
            flux_m1_local[cell.id] += particle.volume * flux_g_tmp * particle.position;
            flux_m2_local[cell.id] +=
                    particle.volume * ((particle.position * particle.position) * flux_g_tmp + flux_h_tmp);
        }
    }
    flux_m0 = mesh.zero_scalar_field();
    flux_m1 = mesh.zero_vector_field();
    flux_m2 = mesh.zero_scalar_field();
    MPI::AllReduce(flux_m0_local, flux_m0);
    MPI::AllReduce(flux_m1_local, flux_m1);
    MPI::AllReduce(flux_m2_local, flux_m2);

    auto m3_local = mesh.zero_vector_field();
    for (auto &cell: mesh.cells) {
        auto dt_v = dt / cell.volume;
        /// tn = n
        auto rho_n = rho_cell[cell.id];
        auto u_n = vel_cell[cell.id];
        auto T_n = T_cell[cell.id];
        auto q_n = q_cell[cell.id];
        auto rhoU_n = rho_n * u_n;
        auto rhoE_n = rho_n * (0.5 * (u_n * u_n) + Cv * T_n);
        auto tau_n = tau_f(rho_n, T_n);
        /// tn = n + 1
        auto rho = rho_n - dt_v * flux_m0[cell.id];
        auto rhoU = rhoU_n - dt_v * flux_m1[cell.id];
        auto rhoE = rhoE_n - dt_v * 0.5 * flux_m2[cell.id];
        auto u = rhoU / rho;
        auto T = (rhoE / rho - (u * u) * 0.5) / Cv;
        auto tau = tau_f(rho, T);
        /// Cell Macro Vars
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = u;
        T_cell[cell.id] = T;
        /// Evolution
        auto cm = half_dt / (half_dt - 2.0 * tau_n);
        auto C = tau / (tau + half_dt);
        auto C_s = half_dt / (2.0 * tau);
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
            auto cq = c * q_n;
            auto g_m = g_maxwell(rho, T, cc);
            auto g_s = g_shakhov(rho, T, cc, cq, g_m);
            auto h_s = h_shakhov(rho, T, cc, cq, g_m);
            /// Evolution
            auto g = C * (g_n + half_dt * (g_s / tau + (g_s_n - g_n) / tau_n) - dt_v * flux_g[p][cell.id]);
            auto h = C * (h_n + half_dt * (h_s / tau + (h_s_n - h_n) / tau_n) - dt_v * flux_h[p][cell.id]);
            /// heat-flux
            m3_local[cell.id] += 0.5 * particle.volume * c * (cc * g + h);
            /// f -> f_bar_plus
            g_cell[p][cell.id] = (1.0 - C_s) * g + C_s * g_s;
            h_cell[p][cell.id] = (1.0 - C_s) * h + C_s * h_s;
        }
        /// Monitor
        if (std::isnan(rho)) {
            logger.warn << "[WARN] cell<" << cell.id + 1 << ">@" << cell.position.str()
                        << " caught rho=nan." << std::endl;
            auto &group = config.get_cell_group(cell, mesh);
            auto mask = group.patch.get_switch("mask");
            if (mask) {
                cell_value_interp(cell.id);
            } else {
                run_state = false;
            }
        }
        if (T < 0.0 or std::isnan(T)) {
            logger.warn << "[WARN] cell<" << cell.id + 1 << ">@" << cell.position.str()
                        << " caught T < 0 or T=nan." << std::endl;
            auto &group = config.get_cell_group(cell, mesh);
            auto mask = group.patch.get_switch("mask");
            if (mask) {
                cell_value_interp(cell.id);
            } else {
                run_state = false;
            }
        }
    }
    MPI::AllReduce(m3_local, q_cell);
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
            List <Scalar> residual_list;
            if (mesh.dimension() == 2) {
                residual_list = {rho_res, T_res, vel_res.x, vel_res.y, q_res.x, q_res.y};
                Utils::print_names_and_values({"Res[Rho]", "Res[T]", "Res[u]", "Res[v]", "Res[qx]", "Res[qy]"},
                                              residual_list);
            } else {
                residual_list = {rho_res, T_res, vel_res.x, vel_res.y, vel_res.z, q_res.x, q_res.y, q_res.z};
                Utils::print_names_and_values(
                        {"Res[Rho]", "Res[T]", "Res[u]", "Res[v]", "Res[w]", "Res[qx]", "Res[qy]", "Res[qz]"},
                        residual_list);
            }
            if (Utils::is_converged(residual_list, residual_limit)) {
                converge_state++;
            } else {
                converge_state = 0;
            }
            if (converge_state >= 10) {
                logger.note << "Step " << step << " - Convergence achieved: residual all below " << residual_limit
                            << std::endl;
                run_state = false;
            }
        }
    }
    MPI::Bcast(run_state);
    MPI_Barrier(MPI_COMM_WORLD);
    /// update Config when step over
    update_config();
    MPI_Barrier(MPI_COMM_WORLD);
}

void CDUGKS_SHAKHOV::cell_value_interp(MESO::ObjectId cell_id) {
    auto rho = fvmMesh::interp_IDW(rho_cell, cell_id);
    auto u = fvmMesh::interp_IDW(vel_cell, cell_id);
    auto T = fvmMesh::interp_IDW(T_cell, cell_id);
    auto q = fvmMesh::interp_IDW(q_cell, cell_id);
    rho_cell[cell_id] = rho;
    vel_cell[cell_id] = u;
    T_cell[cell_id] = T;
    q_cell[cell_id] = q;
    for (int p = 0; p < mpi_task.size; ++p) {
        ObjectId dvs_id = p + mpi_task.start;
        auto &particle = dvs_mesh.cells[dvs_id];
        auto c = particle.position - u;
        auto cc = c * c;
        auto cq = c * q;
        auto gm = g_maxwell(rho, T, cc);
        g_cell[p][cell_id] = g_shakhov(rho, T, cc, cq, gm);
        h_cell[p][cell_id] = h_shakhov(rho, T, cc, cq, gm);
    }
    logger.warn << "[Warn] Mask working on cell<" << cell_id + 1 << "> " << mesh.cells[cell_id].position.str()
                << std::endl;
    Utils::print_names_and_values({"Rho", "u", "v", "w", "T", "qx", "qy", "qz"},
                                  {rho, u.x, u.y, u.z, T, q.x, q.y, q.z});
}

void CDUGKS_SHAKHOV::output() {
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
        return;
    }
    if (step == 0) {
        mesh.output_grid(case_name + "/grid.plt");
    }

    String path, path_latest;
    {
        StringStream ss1, ss2;
        ss1 << case_name << "/step-" << step << "/";
        ss2 << case_name << "/latest/";
        path = ss1.str();
        path_latest = ss2.str();
    }
    Utils::mkdir(path);
    Utils::mkdir(path_latest);

    logger.note << "Output: ";
    auto U = vel_cell.heft(0);
    auto V = vel_cell.heft(1);
    auto W = vel_cell.heft(2);
    auto qx = q_cell.heft(0);
    auto qy = q_cell.heft(1);
    auto qz = q_cell.heft(2);
    if (mesh.dimension() == 2) {
        mesh.output_data(path + "/data.plt",
                         {rho_cell, U, V, T_cell, qx, qy},
                         {"Rho", "U", "V", "T", "qx", "qy"},
                         solution_time);
    } else {
        mesh.output_data(path + "/data.plt",
                         {rho_cell, U, V, W, T_cell, qx, qy, qz},
                         {"Rho", "U", "V", "W", "T", "qx", "qy", "qz"},
                         solution_time);
    }
    logger.note << " tec-files";

    if (output_np) {
        /// output numpy data
        rho_cell.output(path + "Rho");
        vel_cell.output(path + "U");
        T_cell.output(path + "T");
        q_cell.output(path + "q");
        flux_m0.output(path + "fluxM0");
        flux_m1.output(path + "fluxM1");
        flux_m2.output(path + "fluxM2");

        logger.note << ", np-data";
    }
    logger.note << " @ step=" << step << std::endl;
}
