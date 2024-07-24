#include "meso.h"

std::string logo = "\n"
                   "    __  ______________ ____ \n"
                   "   /  |/  / ____/ ___// __ \\\n"
                   "  / /|_/ / __/  \\__ \\/ / / /\n"
                   " / /  / / /___ ___/ / /_/ / \n"
                   "/_/  /_/_____//____/\\____/  \n"
                   "                            \n";

void MESO::help() {
    logger.info << logo << std::endl;
    logger.note << "Usage: \"mpirun [<mpi-params>] ./meso-mpi [<meso-params>]\"" << std::endl;
    logger.info << "   meso-params:\n"
                 "     --case /path/to/case-file       \n"
                 "     --max-step <value:int, default: >          \n"
                 "     --save-interval <value:int, default: 1000>     \n"
                 "     --solver <solver:str>           \n" << std::endl;
}


int MESO::handle_mesh(const std::string &mesh_file, double mesh_scale) {
    auto mesh = MESO::fvmMesh::load_gambit(mesh_file);
    mesh.build_geom(mesh_scale);
    mesh.partition();
    mesh.info();
    mesh.output(mesh_file);

    auto cell_pos = mesh.zero_vector_field();
    cell_pos.MeshCellValueToField([](fvmMesh::Cell &cell) { return cell.position; });
    cell_pos.output(mesh_file + "-cell-position");
    return 0;
}


/// Solver handler
using namespace MESO;
int Solver::solver_state = 0;

#ifdef SOLVER_CDUGKS_H

template<>
void Solver::solver_interrupt<Solver::CDUGKS>(int signum) {
    solver_state = signum;
#pragma omp master
    logger.warn << "Interrupt signal (" << signum << ") received." << std::endl;
}

template<>
int Solver::handle_solver<Solver::CDUGKS>(MESO::ArgParser &parser, MESO::Solver::Config &config) {
    int save_interval = parser.parse_param<int>("save-interval", DefaultValue::save_interval, false);
    MESO::Solver::CDUGKS solver(parser, config);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS>);    // 注册 SIGINT 信号处理
    const int max_step = parser.parse_param("max-step", DefaultValue::max_step, false);
    for (int i = 0; (i < max_step) || (max_step < 0); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI_Barrier(MPI_COMM_WORLD);
            return solver_state;
        }
    }
    logger.note << "Reached the maximum number of iterations: " << max_step << std::endl;
    solver.output();
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#endif

#ifdef SOLVER_CDUGKS_SHAKHOV_H

template<>
void Solver::solver_interrupt<Solver::CDUGKS_SHAKHOV>(int signum) {
    solver_state = signum;
#pragma omp master
    logger.warn << "Interrupt signal (" << signum << ") received." << std::endl;
}

template<>
int Solver::handle_solver<Solver::CDUGKS_SHAKHOV>(MESO::ArgParser &parser, MESO::Solver::Config &config) {
    int save_interval = parser.parse_param<int>("save-interval", DefaultValue::save_interval, false);
    MESO::Solver::CDUGKS_SHAKHOV solver(parser, config);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS_SHAKHOV>);    // 注册 SIGINT 信号处理
    const int max_step = parser.parse_param("max-step", DefaultValue::max_step, false);
    for (int i = 0; (i < max_step) || (max_step < 0); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI_Barrier(MPI_COMM_WORLD);
            return solver_state;
        }
    }
    logger.note << "Reached the maximum number of iterations: " << max_step << std::endl;
    solver.output();
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#endif

