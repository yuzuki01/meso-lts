#ifndef SOLVER_DUGKS_H
#define SOLVER_DUGKS_H

class DUGKS : public BasicSolver {
public:
    fvmMesh::Mesh mesh;
    fvmMesh::Mesh dvs_mesh;
    MPI::MPI_TaskObject mpi_task;

    bool gradient_switch;
    double Ma, Re, CFL;
    double RT, Rho0, L0;
    double tau{}, dt{}, half_dt{};
    int step{};
    double solution_time{};

    /// FieldValue
    Field<Scalar> rho_cell;
    Field<Vector> vel_cell;
    Field<Scalar> rho_cell_res;
    Field<Vector> vel_cell_res;
    DistributionFunction f_cell;
    DistributionFunction f_face;
    DistributionFunction flux_f;

    DUGKS(ArgParser &parser, Config &config);

    void initial();

    Scalar f_maxwell(double rho, const Vector &u, const Vector &particle_velocity) const;

    void reconstruct();

    void fvm_update();

    void do_step();

    void output();
};

#endif //SOLVER_DUGKS_H
