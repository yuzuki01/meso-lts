#include "solver/solver.h"


using namespace MESO::Solver;


CDUGKS::CDUGKS(MESO::ArgParser &parser) : BasicSolver(parser) {
    /// params
    RT = config.get("gas-constant", 0.5) * config.get("ref-temperature", 1.0);
    Rho0 = config.get("ref-density", 1.0);
    L0 = config.get("ref-length", 1.0);
    Ma = config.get("Ma", 0.1);
    Re = config.get("Re", 400.0);
    CFL = config.get("CFL", 0.8);
    gradient_switch = config.get<bool>("gradient-switch", true);
    if (not gradient_switch) {
        logger.warn << "[Warn] gradient-switch was CLOSED." << std::endl;
    }

    mesh = MESO::Mesh::load_gambit(config.get<std::string>("mesh-file", "<mesh-file>"));
    mesh.info();
    config.info();
    dvs_mesh = MESO::Solver::generate_gauss_hermite(mesh.dimension(), 3, RT);
    mpi_task = MPI::DVS_partition(dvs_mesh);
    dvs_mesh.info();

    logger.note << "Solver::CDUGKS loaded." << std::endl;
}

void CDUGKS::initial() {
    step = 0;
    solution_time = 0.0;
    double c_sound = sqrt(RT);
    tau = (L0 * Ma) / (Re * c_sound);
    dt = CFL * ((mesh.min_cell_size / 2.0) / dvs_mesh.max_cell_magnitude);
    half_dt = 0.5 * dt;

    m0_cell = Field<Scalar>(mesh, cell_field_flag);
    m0_cell_n = Field<Scalar>(mesh, cell_field_flag);
    m0_cell_res = Field<Scalar>(mesh, cell_field_flag);
    m1_cell = Field<Vector>(mesh, cell_field_flag);
    m1_cell_n = Field<Vector>(mesh, cell_field_flag);
    m1_cell_res = Field<Vector>(mesh, cell_field_flag);
    m0_face = Field<Scalar>(mesh, face_field_flag);
    m1_face = Field<Vector>(mesh, face_field_flag);
    f_cell.resize(mpi_task.size, Field<Scalar>(mesh, cell_field_flag));
    f_face.resize(mpi_task.size, Field<Scalar>(mesh, face_field_flag));

#pragma omp parallel for default(none)
    for (auto &cell: mesh.cells) {
        auto &group = config.get_cell_group(cell, mesh);
        double m0 = 0.0;
        Vector m1(0.0, 0.0, 0.0);
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            double f = f_maxwell(group.density, group.velocity, particle.position);
            f_cell[p][cell.id] = f;
            m0 += particle.volume * f;
            m1 += particle.volume * f * particle.position;
        }
        m0_cell_n[cell.id] = m0;
        m1_cell_n[cell.id] = m1;
    }
    /// X_cell_n as X_cell_local here
    MPI::ReduceAll(m0_cell_n, m0_cell);
    MPI::ReduceAll(m1_cell_n, m1_cell);
    MPI_Barrier(MPI_COMM_WORLD);
    m0_cell_n.values = m0_cell.values;
    m1_cell_n.values = m1_cell.values;
    is_crashed = false;
    logger.note << "Initialization finished." << std::endl;
}

MESO::Scalar CDUGKS::f_maxwell(double m0, const Vector &m1, const Vector &particle_velocity) const {
    double over_RT = 1.0 / RT;
    Vector flow_velocity = m1 / m0;
    double ku_RT = (flow_velocity * particle_velocity) * over_RT;
    double uu_RT = (flow_velocity * flow_velocity) * over_RT;
    return m0 * (1.0 + ku_RT + (ku_RT * ku_RT - uu_RT) * 0.5);
}

void CDUGKS::reconstruct() {
    /// get f_bar on face
    for (int p = 0; p < mpi_task.size; ++p) {
        ObjectId dvs_id = p + mpi_task.start;
        auto &particle = dvs_mesh.cells[dvs_id];
        /// cell gradient
        Field <Vector> grad_f = f_cell[p].gradient(gradient_switch);  // Field<Scalar>.gradient 内置 openMP 并行
        /// interp to face
#pragma omp parallel for shared(particle, grad_f, p) default(none)
        for (auto &face: mesh.faces) {
            Vector &nv = face.normal_vector[0];
            auto &cell = (nv * particle.position >= 0.0) ? mesh.cells[face.cell_id[0]] : mesh.cells[face.cell_id[1]];
            f_face[p][face.id] = f_cell[p][cell.id]
                                 + (face.position - cell.position - particle.position * half_dt) * grad_f[cell.id];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    {
        Field <Scalar> m0_local(mesh, face_field_flag);
        Field <Vector> m1_local(mesh, face_field_flag);
#pragma omp parallel for shared(m0_local, m1_local) default(none)
        for (auto &face: mesh.faces) {
            double m0 = 0.0;
            Vector m1(0.0, 0.0, 0.0);
            /// face macro vars
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                double f = f_face[p][face.id];
                m0 += particle.volume * f;
                m1 += particle.volume * f * particle.position;
            }
            m0_local[face.id] = m0;
            m1_local[face.id] = m1;
        }
        MPI::ReduceAll(m0_local, m0_face);
        MPI::ReduceAll(m1_local, m1_face);
    }
    MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel for default(none)
    for (auto &face: mesh.faces) {
        /// get original f on face
        double c_eq = half_dt / (2.0 * tau + half_dt);
        for (int p = 0; p < mpi_task.size; ++p) {
            // openMP parallel runs on mesh
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            double f_eq = f_maxwell(m0_face[face.id], m1_face[face.id], particle.position);
            f_face[p][face.id] = (1.0 - c_eq) * f_face[p][face.id] + c_eq * f_eq;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    {
        Field <Scalar> rho_w_local(mesh, face_field_flag);
        Field <Scalar> rho_w0_local(mesh, face_field_flag);
        Field <Scalar> rho_w_global(mesh, face_field_flag);
        Field <Scalar> rho_w0_global(mesh, face_field_flag);

#pragma omp parallel for shared(rho_w_local, rho_w0_local) default(none)
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
                            f_face[p][face.id] = f_maxwell(mark.density, mark.density * mark.velocity,
                                                           particle.position);
                        }
                    }
                    break;
                case outlet:
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        if (particle.position * nv >= 0.0) {
                            f_face[p][face.id] = f_maxwell(m0_cell[neighbor.id], m1_cell[neighbor.id],
                                                           particle.position);
                        }
                    }
                    break;
                case wall: {
                    double rho_w, rho_w0;
                    rho_w = rho_w0 = 0.0;
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        double kn = particle.position * nv;
                        if (kn >= 0.0) {
                            double f_eq = f_maxwell(1.0, mark.velocity, particle.position);
                            rho_w0 += kn * particle.volume * f_eq;
                            f_face[p][face.id] = f_eq;
                        } else {
                            rho_w -= kn * particle.volume * f_face[p][face.id];
                        }
                    }
                    rho_w_local[face.id] = rho_w;
                    rho_w0_local[face.id] = rho_w0;
                }
                case fluid_interior:
                default:
                    break;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI::ReduceAll(rho_w_local, rho_w_global);
        MPI::ReduceAll(rho_w0_local, rho_w0_global);
        MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel for shared(rho_w_global, rho_w0_global) default(none)
        for (auto &face: mesh.faces) {
            auto &mark = config.get_face_group(face, mesh);
            switch (mark.type) {
                case wall: {
                    auto &nv = face.normal_vector[1];
                    for (int p = 0; p < mpi_task.size; ++p) {
                        ObjectId dvs_id = p + mpi_task.start;
                        auto &particle = dvs_mesh.cells[dvs_id];
                        double kn = particle.position * nv;
                        if (kn >= 0.0) {
                            f_face[p][face.id] *= rho_w_global[face.id] / rho_w0_global[face.id];
                        }
                    }
                }
                    break;
                default:
                    break;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}


void CDUGKS::fvm_update() {
    {
        Field <Scalar> flux_m0_local(mesh, cell_field_flag);
        Field <Vector> flux_m1_local(mesh, cell_field_flag);
        Field <Scalar> flux_m0_global(mesh, cell_field_flag);
        Field <Vector> flux_m1_global(mesh, cell_field_flag);

#pragma omp parallel for shared(flux_m0_local, flux_m1_local) default(none)
        for (auto &cell: mesh.cells) {
            double flux_m0 = 0.0;
            Vector flux_m1(0.0, 0.0, 0.0);
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                double flux = 0.0;
                for (auto face_id: cell.face_id) {
                    auto &face = mesh.faces[face_id];
                    auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                    flux += (particle.position * nv) * face.area * f_face[p][face.id];
                }
                flux_m0 += particle.volume * flux;
                flux_m1 += particle.volume * flux * particle.position;
            }
            flux_m0_local[cell.id] = flux_m0;
            flux_m1_local[cell.id] = flux_m1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI::ReduceAll(flux_m0_local, flux_m0_global);
        MPI::ReduceAll(flux_m1_local, flux_m1_global);
        MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel for shared(flux_m0_global, flux_m1_global) default(none)
        for (auto &cell: mesh.cells) {
            double dt_v = dt / cell.volume;
            m0_cell[cell.id] = m0_cell_n[cell.id] - dt_v * flux_m0_global[cell.id];
            m1_cell[cell.id] = m1_cell_n[cell.id] - dt_v * flux_m1_global[cell.id];
            double C = tau / (tau + half_dt);
            double cm = half_dt / (half_dt - 2.0 * tau);
            double c_eq = half_dt / (2.0 * tau);
            for (int p = 0; p < mpi_task.size; ++p) {
                ObjectId dvs_id = p + mpi_task.start;
                auto &particle = dvs_mesh.cells[dvs_id];
                double flux = 0.0;
                for (auto face_id: cell.face_id) {
                    auto &face = mesh.faces[face_id];
                    auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                    flux += (particle.position * nv) * face.area * f_face[p][face.id];
                }
                double f_eq_n = f_maxwell(m0_cell_n[cell.id], m1_cell_n[cell.id], particle.position);
                double f_n = cm * f_eq_n + (1.0 - cm) * f_cell[p][cell.id];
                double f_eq = f_maxwell(m0_cell[cell.id], m1_cell[cell.id], particle.position);
                // 直接演化得 n + 1 时刻的 f_bp
                f_cell[p][cell.id] = (1.0 - c_eq) *
                                     (C * (f_n + half_dt * (f_eq / tau + (f_eq_n - f_n) / tau) - dt_v * flux)) +
                                     c_eq * f_eq;
            }
            m0_cell_n[cell.id] = m0_cell[cell.id];
            m1_cell_n[cell.id] = m1_cell[cell.id];
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void CDUGKS::do_step() {
    reconstruct();
    fvm_update();
    step++;
    solution_time += dt;
    if (MPI::rank == MPI::main_rank) {
        if (step % residual_interval == 0) {
            double m0_res = residual(m0_cell_res, m0_cell);
            Vector m1_res = residual(m1_cell_res, m1_cell);
            logger.note << "step: " << step << std::endl;
            if (mesh.dimension() == 2) {
                Utils::print_names_and_values({"Res[m0]", "Res[m1x]", "Res[m1y]"},
                                              {m0_res, m1_res.x, m1_res.y});
            } else {
                Utils::print_names_and_values({"Res[m0]", "Res[m1x]", "Res[m1y]", "Res[m1z]"},
                                              {m0_res, m1_res.x, m1_res.y, m1_res.z});
            }
        }
    }
}

void CDUGKS::output() {
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
        return;
    }
    std::stringstream file_name;
    file_name << "./" << case_name << "/result-" << step;

    auto m1x = m1_cell.heft(0);
    auto m1y = m1_cell.heft(1);
    auto m1z = m1_cell.heft(2);
    mesh.output(file_name.str(),
                {"m0", "m1x", "m1y", "m1z"},
                {&m0_cell, &m1x, &m1y, &m1z}, step, solution_time);
}
