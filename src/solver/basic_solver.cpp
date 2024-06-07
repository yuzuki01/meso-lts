#include "solver/solver.h"

using namespace MESO;

Solver::BasicSolver::BasicSolver(ArgParser &parser) :
        parser(parser), config(parser.parse_param<std::string>("case", "<case-file>", true)) {
    case_name = config.get<std::string>("case-name", "unnamed", false);
    output_np = config.get<bool>("output-np", false, false);
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
    }
}


/// Solver handler
int Solver::solver_state = 0;

#ifdef SOLVER_CDUGKS_H
template<>
void Solver::solver_interrupt<Solver::CDUGKS>(int signum) {
    solver_state = signum;
#pragma omp master
    logger.warn << "Interrupt signal (" << signum << ") received." << std::endl;
}

template<>
int Solver::handle_solver<Solver::CDUGKS>(int *p_argc, char ***p_argv, MESO::ArgParser &parser) {
    /// MPI init
    MESO::MPI::Initialize(p_argc, p_argv, parser.parse_param<int>("parallel", 1, true));
    int save_interval = parser.parse_param<int>("save-interval", 1000, false);
    MESO::Solver::CDUGKS solver(parser);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MESO::MPI::Finalize();
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS>);    // 注册 SIGINT 信号处理
    for (int i = 0; i < parser.parse_param("max-step", 10000, false); ++i) {
        solver.do_step();
        if (solver.is_crashed) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI::Finalize();
            return solver_state;
        }
    }
    solver.output();
    MPI::Finalize();
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
int Solver::handle_solver<Solver::CDUGKS_SHAKHOV>(int *p_argc, char ***p_argv, MESO::ArgParser &parser) {
    /// MPI init
    MESO::MPI::Initialize(p_argc, p_argv, parser.parse_param<int>("parallel", 1, true));
    int save_interval = parser.parse_param<int>("save-interval", 1000, false);
    MESO::Solver::CDUGKS_SHAKHOV solver(parser);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        MESO::MPI::Finalize();
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS_SHAKHOV>);    // 注册 SIGINT 信号处理
    for (int i = 0; i < parser.parse_param("max-step", 10000, false); ++i) {
        solver.do_step();
        if (solver.is_crashed) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            MPI::Finalize();
            return solver_state;
        }
    }
    solver.output();
    MPI::Finalize();
    return 0;
}
#endif
