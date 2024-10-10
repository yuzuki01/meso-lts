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
    logger.info << "\n\n./meso-mpi [Options]\n"
                   "    --mesh <mesh-file>";
}


int MESO::handle_mesh(const std::string &mesh_file, double mesh_scale) {
    auto mesh = MESO::fvmMesh::load_gambit(mesh_file);
    mesh.build_geom(mesh_scale);
#ifdef _METIS_H_
    mesh.partition();
#endif
    mesh.info();
    mesh.output(mesh_file);
    return 0;
}

/// Solver handler
using namespace MESO;
int Solver::solver_state = 0;

#ifdef SOLVER_DUGKS_H

template<>
int Solver::handle_solver<Solver::DUGKS>(MESO::ArgParser &parser, MESO::Solver::Config &config) {
    MESO::Solver::DUGKS solver(parser, config);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }

    for (int i = 0; (i < solver.max_step()) || (solver.max_step() < 0); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % solver.write_interval() == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI_Barrier(MPI_COMM_WORLD);
            return solver_state;
        }
    }
    logger.note << "Reached stop step: " << solver.max_step() << std::endl;
    solver.output();
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#endif

#ifdef SOLVER_CDUGKS_H

template<>
int Solver::handle_solver<Solver::CDUGKS>(MESO::ArgParser &parser, MESO::Solver::Config &config) {
    MESO::Solver::CDUGKS solver(parser, config);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }
    
    for (int i = 0; (i < solver.max_step()) || (solver.max_step() < 0); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % solver.write_interval() == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI_Barrier(MPI_COMM_WORLD);
            return solver_state;
        }
    }
    logger.note << "Reached stop step: " << solver.max_step() << std::endl;
    solver.output();
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#endif

#ifdef SOLVER_CDUGKS_SHAKHOV_H

template<>
int Solver::handle_solver<Solver::CDUGKS_SHAKHOV>(MESO::ArgParser &parser, MESO::Solver::Config &config) {
    MESO::Solver::CDUGKS_SHAKHOV solver(parser, config);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }

    for (int i = 0; (i < solver.max_step()) || (solver.max_step() < 0); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % solver.write_interval() == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI_Barrier(MPI_COMM_WORLD);
            return solver_state;
        }
    }
    logger.note << "Reached stop step: " << solver.max_step() << std::endl;
    solver.output();
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#endif
