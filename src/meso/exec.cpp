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


void MESO::version() {
    // 1. 编译日期和时间
    logger.info << "Compile Date/Time:    " << COMPILE_DATE << " " << COMPILE_TIME << std::endl;

    // 2. 编译器信息
    logger.info << "Compiler:             " << COMPILER_NAME << std::endl;
    logger.info << "Compiler Version:     " << COMPILER_VERSION_STR << " (Numeric: " << COMPILER_VERSION << ")" << std::endl;

    // 3. 系统架构（可选扩展）
#if defined(_WIN32) || defined(_WIN64)
    #if defined(_WIN64)
        logger.info << "System Arch:          x86_64 (Windows 64-bit)" << std::endl;
    #else
        logger.info << "System Arch:          x86 (Windows 32-bit)" << std::endl;
    #endif
#elif defined(__x86_64__) || defined(_M_X64)
    logger.info << "System Arch:          x86_64 (Linux/macOS 64-bit)" << std::endl;
#elif defined(__i386__) || defined(_M_IX86)
    logger.info << "System Arch:          x86 (Linux/macOS 32-bit)" << std::endl;
#elif defined(__arm64__) || defined(_M_ARM64)
    logger.info << "System Arch:          ARM64" << std::endl;
#else
    logger.info << "System Arch:          Unknown" << std::endl;
#endif

    // 4. C++ 标准版本（可选扩展）
#if __cplusplus == 201703L
    logger.info << "C++ Standard:         C++17" << std::endl;
#elif __cplusplus == 201402L
    logger.info << "C++ Standard:         C++14" << std::endl;
#elif __cplusplus == 201103L
    logger.info << "C++ Standard:         C++11" << std::endl;
#elif __cplusplus == 199711L
    logger.info << "C++ Standard:         C++98" << std::endl;
#else
    logger.info << "C++ Standard:         C++20+ (Version: " << __cplusplus << ")" << std::endl;
#endif
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
