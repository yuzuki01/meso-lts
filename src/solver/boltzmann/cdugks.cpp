#include "solver/solver.h"


using namespace MESO::Solver;


CDUGKS::CDUGKS(MESO::ArgParser &parser, Config &config) : BasicSolver(parser, config) {
    /// params
    RT = config.get("gas-constant", 0.5) * config.get("ref-temperature", 1.0);
    Rho0 = config.get("ref-density", 1.0);
    L0 = config.get("ref-length", 1.0);
    Ma = config.get("Ma", 0.1);
    Re = config.get("Re", 400.0);
    dt = config.get("timeStep", 2.5e-3);
    gradient_switch = config.get<bool>("gradient-switch", true);
    if (not gradient_switch) {
        logger.warn << "[Warn] gradient-switch was CLOSED." << std::endl;
    }

    mesh = MESO::fvmMesh::load_gambit(config.get<String>("mesh-file", "<mesh-file>"),
                                      mesh_scale);
    mesh.info();
    config.info();

    auto dvs_file = config.get<String>("dvs-file", "<mesh-file>");
    auto dvs_type = config.get<String>("dvs-type", "<dvs-type>");
    if (dvs_type == "Gauss-Hermite") {
        DVS::GaussHermiteParams dvs_params(mesh.dimension(), RT);
        dvs_mesh = MESO::DVS::generate_dvs(dvs_file,
                                           dvs_params);
    } else {
        logger.error << "dvs-mesh type error." << std::endl;
        throw std::invalid_argument("dvs-mesh type error.");
    }
    mpi_task = MPI::DVS_partition(dvs_mesh);
    dvs_mesh.info();

    logger.note << "Solver::CDUGKS loaded." << std::endl;
}

void CDUGKS::initial() {
    step = 0;
    converge_state = 0;
    solution_time = 0.0;
    Scalar c_sound = sqrt(RT);
    tau = (L0 * Ma) / (Re * c_sound);
    half_dt = 0.5 * dt;

    rho_cell = mesh.zero_scalar_field();
    vel_cell = mesh.zero_vector_field();
    f_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    f_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));
    /// flux
    flux_f.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));

    if (config.get("read-np-data", false, false)) {
        read_np_data<Scalar>(case_name + "/latest/Rho", rho_cell);
        read_np_data<Vector>(case_name + "/latest/U", vel_cell);
        for (auto &cell: mesh.cells) {
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                f_cell[p][cell.id] = f_maxwell(rho_cell[cell.id], vel_cell[cell.id], particle.position);
            }
        }
    } else {
        auto m0_local = mesh.zero_scalar_field();
        auto m1_local = mesh.zero_vector_field();
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
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                Scalar f = f_maxwell(rho_patch, u_patch, particle.position);
                f_cell[p][cell.id] = f;
                m0_local[cell.id] += particle.volume * f;
                m1_local[cell.id] += particle.volume * f * particle.position;
            }
        }
        auto m0 = mesh.zero_scalar_field();
        auto m1 = mesh.zero_vector_field();
        MPI::AllReduce(m0_local, m0);
        MPI::AllReduce(m1_local, m1);
        for (auto &cell: mesh.cells) {
            auto rho = m0[cell.id];
            auto rhoU = m1[cell.id];
            auto u = rhoU / rho;
            rho_cell[cell.id] = rho;
            vel_cell[cell.id] = u;
        }
    }

    /// residual
    rho_cell_res = rho_cell;
    vel_cell_res = vel_cell;

    run_state = true;
    logger.note << "Initialization finished." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

MESO::Scalar CDUGKS::f_maxwell(Scalar rho, const Vector &u, const Vector &particle_velocity) const {
    Scalar over_RT = 1.0 / RT;
    Scalar ku_RT = (u * particle_velocity) * over_RT;
    Scalar uu_RT = (u * u) * over_RT;
    return rho * (1.0 + ku_RT + (ku_RT * ku_RT - uu_RT) * 0.5);
}

void CDUGKS::reconstruct() {
    /// get f_bar on face
    for (int p = 0; p < mpi_task.size; ++p) {
        ObjectId dvs_id = p + mpi_task.start;
        auto &particle = dvs_mesh.cells[dvs_id];
        /// cell gradient
        Field <Vector> grad_f = mesh.zero_vector_field(cell_field_flag);
        if (gradient_switch) grad_f = fvmMesh::grad(f_cell[p]);
        /// interp to face
        for (auto &face: mesh.faces) {
            Vector &nv = face.normal_vector[0];
            auto &own = mesh.cells[face.cell_id[0]];
            auto &nei = mesh.cells[face.cell_id[1]];
            auto nv_xi = nv * particle.position;
            if (nv_xi >= VSMALL) {
                // from own
                f_face[p][face.id] = f_cell[p][own.id]
                                     + (face.position - own.position - particle.position * half_dt) * grad_f[own.id];
            } else if (nv_xi < -VSMALL) {
                // from nei
                f_face[p][face.id] = f_cell[p][nei.id]
                                     + (face.position - nei.position - particle.position * half_dt) * grad_f[nei.id];
            } else {
                f_face[p][face.id] = 0.5 * (f_cell[p][own.id] + f_cell[p][nei.id]
                                            + (face.position - own.position - particle.position * half_dt) *
                                              grad_f[own.id]
                                            + (face.position - nei.position - particle.position * half_dt) *
                                              grad_f[nei.id]);
            }
        }
    }
    {
        auto m0_local = mesh.zero_scalar_field(face_field_flag);
        auto m1_local = mesh.zero_vector_field(face_field_flag);
        /// face macro vars
        for (auto &face: mesh.faces) {
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                Scalar f = f_face[p][face.id];
                m0_local[face.id] += particle.volume * f;
                m1_local[face.id] += particle.volume * f * particle.position;
            }
        }
        auto m0 = mesh.zero_scalar_field(face_field_flag);
        auto m1 = mesh.zero_vector_field(face_field_flag);
        MPI::AllReduce(m0_local, m0);
        MPI::AllReduce(m1_local, m1);
        Scalar c_eq = half_dt / (2.0 * tau + half_dt);
        /// get original f on face
        for (auto &face: mesh.faces) {
            auto rho = m0[face.id];
            auto rhoU = m1[face.id];
            auto u = rhoU / rho;
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                auto f_eq = f_maxwell(rho, u, particle.position);
                f_face[p][face.id] = (1.0 - c_eq) * f_face[p][face.id] + c_eq * f_eq;
            }
        }
    }
    {
        Map <ObjectId> wall_rho_map;
        List <Scalar> wall_rho_local, wall_rho0_local;

        /// boundary
        for (auto &face: mesh.faces) {
            auto &mark = config.get_face_group(face, mesh);
            auto &nv = face.normal_vector[1];
            auto &neighbor = mesh.cells[face.cell_id[0]];
            switch (mark.type) {
                case BoundaryType::farfield_inlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            f_face[p][face.id] = f_maxwell(mark.patch.get_scalar("density"),
                                                           mark.patch.get_vector("velocity"),
                                                           particle.position);
                        }
                    }
                    break;
                case BoundaryType::farfield_outlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= VSMALL) {
                            f_face[p][face.id] = f_maxwell(rho_cell[neighbor.id], vel_cell[neighbor.id],
                                                           particle.position);
                        }
                    }
                    break;
                case BoundaryType::wall: {
                    Scalar rho_w, rho_w0;
                    rho_w = rho_w0 = 0.0;
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        Scalar kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            Scalar f_eq = f_maxwell(1.0, mark.patch.get_vector("velocity"),
                                                    particle.position);
                            rho_w0 += kn * particle.volume * f_eq;
                        } else {
                            rho_w -= kn * particle.volume * f_face[p][face.id];
                        }
                    }
                    wall_rho_map[face.id] = int(wall_rho_local.size());
                    wall_rho_local.push_back(rho_w);
                    wall_rho0_local.push_back(rho_w0);
                }
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
                case wall: {
                    auto &nv = face.normal_vector[1];
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        Scalar kn = particle.position * nv;
                        if (kn >= VSMALL) {
                            int wall_rho_list_id = wall_rho_map[face.id];
                            auto rho_w = wall_rho_global[wall_rho_list_id] / wall_rho0_global[wall_rho_list_id];
                            f_face[p][face.id] = f_maxwell(rho_w, mark.patch.get_vector("velocity"),
                                                           particle.position);
                        } else {
                            if (kn < -VSMALL) continue;
                            int wall_rho_list_id = wall_rho_map[face.id];
                            auto rho_w = wall_rho_global[wall_rho_list_id] / wall_rho0_global[wall_rho_list_id];
                            f_face[p][face.id] = 0.5 * (
                                    f_maxwell(rho_w,
                                              mark.patch.get_vector("velocity"),
                                              particle.position)
                                    + f_face[p][face.id]
                            );
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


void CDUGKS::fvm_update() {
    /// Flux
    auto flux_m0_local = mesh.zero_scalar_field();
    auto flux_m1_local = mesh.zero_vector_field();

    for (auto &cell: mesh.cells) {
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            Scalar flux = 0.0;
            for (auto face_id: cell.face_id) {
                auto &face = mesh.faces[face_id];
                auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                flux += (particle.position * nv) * face.area * f_face[p][face.id];
            }
            flux_f[p][cell.id] = flux;
            flux_m0_local[cell.id] += particle.volume * flux;
            flux_m1_local[cell.id] += particle.volume * flux * particle.position;
        }
    }
    auto flux_m0 = mesh.zero_scalar_field();
    auto flux_m1 = mesh.zero_vector_field();
    MPI::AllReduce(flux_m0_local, flux_m0);
    MPI::AllReduce(flux_m1_local, flux_m1);

    auto cm = half_dt / (half_dt - 2.0 * tau);
    auto C = tau / (tau + half_dt);
    auto c_eq = half_dt / (2.0 * tau);

    for (auto &cell: mesh.cells) {
        auto dt_v = dt / cell.volume;
        /// tn = n
        auto rho_n = rho_cell[cell.id];
        auto u_n = vel_cell[cell.id];
        /// tn = n + 1
        auto rho = rho_n - dt_v * flux_m0[cell.id];
        auto u = (rho_n * u_n - dt_v * flux_m1[cell.id]) / rho_n;
        /// Cell Macro Vars
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = u;
        /// Evolution
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            auto f_eq_n = f_maxwell(rho_n, u_n, particle.position);
            auto f_n = cm * f_eq_n + (1.0 - cm) * f_cell[p][cell.id];
            auto f_eq = f_maxwell(rho, u, particle.position);
            // 直接演化得 n + 1 时刻的 f_bp
            f_cell[p][cell.id] = (1.0 - c_eq) *
                                 (C * (f_n + half_dt * (f_eq / tau + (f_eq_n - f_n) / tau) -
                                       dt_v * flux_f[p][cell.id])) +
                                 c_eq * f_eq;
        }
    }
}

void CDUGKS::do_step() {
    reconstruct();
    fvm_update();
    step++;
    solution_time += dt;
    if (MPI::rank == MPI::main_rank) {
        if (step % residual_interval == 0) {
            Scalar m0_res = residual(rho_cell_res, rho_cell);
            Vector m1_res = residual(vel_cell_res, vel_cell);
            logger.note << "step: " << step << std::endl;
            List <Scalar> residual_list;
            if (mesh.dimension() == 2) {
                residual_list = {m0_res, m1_res.x, m1_res.y};
                Utils::print_names_and_values({"Res[Rho]", "Res[U]", "Res[V]"},
                                              residual_list);
            } else {
                residual_list = {m0_res, m1_res.x, m1_res.y, m1_res.z};
                Utils::print_names_and_values({"Res[Rho]", "Res[U]", "Res[V]", "Res[W]"},
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
}

void CDUGKS::output() {
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
    if (mesh.dimension() == 2) {
        mesh.output_data(path + "/data.plt",
                         {rho_cell, U, V},
                         {"Rho", "U", "V"},
                         solution_time);
    } else {
        mesh.output_data(path + "/data.plt",
                         {rho_cell, U, V, W},
                         {"Rho", "U", "V", "W"},
                         solution_time);
    }
    logger.note << " tec-files";

    if (output_np) {
        /// output numpy data
        rho_cell.output(path + "Rho");
        vel_cell.output(path + "U");

        logger.note << ", np-data";
    }
    {
        rho_cell.output(path_latest + "Rho");
        vel_cell.output(path_latest + "U");
    }
    logger.note << " @ step=" << step << std::endl;
}
