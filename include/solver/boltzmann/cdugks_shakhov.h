#ifndef SOLVER_CDUGKS_SHAKHOV_H
#define SOLVER_CDUGKS_SHAKHOV_H

class CDUGKS_SHAKHOV : public BasicSolver {
public:
    fvmMesh::Mesh mesh;
    fvmMesh::Mesh dvs_mesh;
    MPI::MPI_TaskObject mpi_task;

    bool gradient_switch, limiter_switch;
    double Ma, Re{}, Kn, Pr, CFL;
    int D, K;
    double R, Rho0, T0, L0, mfp{};
    double gamma{}, Cv{}, vhs_index, venkata_k;
    double miu0{}, dt{}, half_dt{};
    int step{};
    double solution_time{};

    /// FieldValue
    Field<Scalar> rho_cell, T_cell;
    Field<Vector> vel_cell, q_cell;

    Field<Scalar> rho_cell_res, T_cell_res;
    Field<Vector> vel_cell_res, q_cell_res;

    DistributionFunction g_cell, h_cell;
    DistributionFunction g_face, h_face;
    DistributionFunction flux_g, flux_h;

    CDUGKS_SHAKHOV(ArgParser &parser, Config &config);

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

    template<class FieldDataType>
    void read_np_data(const std::string &file, Field<FieldDataType> &field);
};

#endif //SOLVER_CDUGKS_SHAKHOV_H
