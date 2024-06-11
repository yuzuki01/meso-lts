#include "solver/solver.h"

using namespace MESO;

Solver::BasicSolver::BasicSolver(ArgParser &parser) :
        parser(parser), config(parser.parse_param<std::string>("case", "<case-file>", true)) {
    case_name = config.get<std::string>("case-name", "unnamed", false);
    output_np = config.get<bool>("output-np", false, false);
    if (Utils::mkdir(case_name) != 0) {
        logger.warn << "Solver::output() cannot mkdir: " << case_name << std::endl;
    }
    residual_limit = config.get<double>("residual-limit", 1e-6, false);
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
    int save_interval = parser.parse_param<int>("save-interval", 1000, false);
    MESO::Solver::CDUGKS solver(parser);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS>);    // 注册 SIGINT 信号处理
    for (int i = 0; i < parser.parse_param("max-step", 10000, false); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            return solver_state;
        }
    }
    solver.output();
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
    int save_interval = parser.parse_param<int>("save-interval", 1000, false);
    MESO::Solver::CDUGKS_SHAKHOV solver(parser);
    solver.initial();
    solver.output();
    if (parser.parse_switch("not-run")) {
        return 0;
    }
    std::signal(SIGINT, Solver::solver_interrupt<Solver::CDUGKS_SHAKHOV>);    // 注册 SIGINT 信号处理
    for (int i = 0; i < parser.parse_param("max-step", 10000, false); ++i) {
        solver.do_step();
        if (not solver.get_run_state()) break;
        if (solver.step % save_interval == 0) solver.output();
        if (solver_state != 0) {
            solver.output();
            return solver_state;
        }
    }
    solver.output();
    return 0;
}

#endif
