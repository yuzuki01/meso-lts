#ifndef SOLVER_CDUGKS_H
#define SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
public:
    Mesh::Mesh mesh;
    Mesh::Mesh dvs_mesh;
    MPI::MPI_TaskObject mpi_task;

    bool is_crashed{}, gradient_switch;
    double Ma, Re, CFL;
    double RT, Rho0, L0;
    double tau{}, dt{}, half_dt{};
    int step{};
    double solution_time{};

    /// FieldValue
    Field<Scalar> m0_cell, m0_cell_n, m0_face;
    Field<Vector> m1_cell, m1_cell_n, m1_face;
    DistributionFunction f_cell;
    DistributionFunction f_face;

    explicit CDUGKS(ArgParser &parser);

    void initial();

    Scalar f_maxwell(double m0, const Vector &m1, const Vector &particle_velocity) const;

    void reconstruct();

    void fvm_update();

    void do_step();

    void output();
};

#endif //SOLVER_CDUGKS_H
