#ifndef SOLVER_CDUGKS_SHAKHOV_H
#define SOLVER_CDUGKS_SHAKHOV_H

class CDUGKS_SHAKHOV : public BasicSolver {
public:
    Mesh::Mesh mesh;
    Mesh::Mesh dvs_mesh;
    MPI::MPI_TaskObject mpi_task;

    bool is_crashed{}, gradient_switch;
    double Ma, Re{}, Kn, Pr, CFL;
    int D, K;
    double R, Rho0, T0, L0;
    double gamma{}, Cv{}, vhs_omega, vhs_index;
    double miu0{}, dt{}, half_dt{};
    int step{};
    double solution_time{};

    /// FieldValue
    Field<Scalar> rho_cell, T_cell, tau_cell,
            rho_cell_n, T_cell_n, tau_cell_n;
    Field<Vector> vel_cell, q_cell,
            vel_cell_n, q_cell_n;

    Field<Scalar> rho_cell_res, T_cell_res; // res
    Field<Vector> vel_cell_res, q_cell_res; // res

    Field<Scalar> rho_face, T_face, tau_face;
    Field<Vector> vel_face, q_face;

    DistributionFunction g_cell, h_cell;
    DistributionFunction g_face, h_face;
    DistributionFunction flux_g, flux_h;

    explicit CDUGKS_SHAKHOV(ArgParser &parser);

    void initial();

    inline Scalar tau_f(double rho, double t) const;

    inline Scalar g_maxwell(double rho, double t, double cc) const;

    inline Scalar g_shakhov(double rho, double t, double cc, double cq, double gm) const;

    inline Scalar h_maxwell(double t, double gm) const;

    inline Scalar h_shakhov(double rho, double t, double cc, double cq, double gm) const;

    void reconstruct();

    void fvm_update();

    void do_step();

    void output();
};

#endif //SOLVER_CDUGKS_SHAKHOV_H
