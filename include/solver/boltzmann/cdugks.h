#ifndef SOLVER_CDUGKS_H
#define SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
public:
    Mesh::Mesh mesh;
    Mesh::Mesh dvs_mesh;
    /// mpi
    MPI_Task_List dvs_partition; // [node_rank]{start_index, num}

    bool is_crashed{}, gradient_switch;
    double Ma, Re, CFL;
    double RT, Rho0, L0;
    double tau{}, dt{}, half_dt{};
    int step{};
    double solution_time{};

    /// FieldValue
    Field<Scalar> rho_cell;
    Field<Vector> vel_cell;
    DistributionFunction f_cell;
    DistributionFunction f_face;

    explicit CDUGKS(ArgParser &parser);

    void initial();

    Scalar f_maxwell(double density, const Vector &flow_velocity, const Vector &particle_velocity) const;

    void reconstruct();

    void fvm_update();

    void do_step();

    void output();
};

#endif //SOLVER_CDUGKS_H
