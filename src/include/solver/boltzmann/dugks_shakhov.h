#ifndef SOLVER_DUGKS_SHAKHOV
#define SOLVER_DUGKS_SHAKHOV

class DUGKS_SHAKHOV : public BasicSolver {
public:
    /// Constructor
    using Scheme = DUGKS_SHAKHOV;
    explicit DUGKS_SHAKHOV(ConfigReader &_config, ArgParser &_parser);
    /// Check Point
    CheckPoint<DUGKS_SHAKHOV> check_point;
    /// Mesh
    MESH::ListMesh phy_mesh{MESH_TYPE_NORMAL, config.phy_mesh};
    MESH::ListMesh dvs_mesh{MESH_TYPE_NO_FACE, config.dvs_mesh};
    /// Physical
    double Re{}, Ma, Kn, Pr;
    double R, T0, Rho0, L0, vhs_index, CFL;
    double dt{}, half_dt{}, miu0{}, Cv{}, gamma{};
    int D, K;
    bool limiter_switch, zero_gradient, boundary_zero_gradient;
    double ventaka_k{};

    /// Physical Formula
    inline double tau_f(double density, double temperature) const;
    inline double g_maxwell(double density, double temperature, const Vec3D &c) const;
    inline double g_shakhov(double density, double temperature, const Vec3D &c, const Vec3D &heat_flux) const;
    inline double h_maxwell(double g_m, double temperature) const;
    inline double h_shakhov(double density, double temperature, const Vec3D &c, const Vec3D &heat_flux) const;

    using DistributionFunction = std::vector<double>;
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> &mesh_cell;
        /// 最小二乘法
        LeastSquare lsp;
        /// Limiter
        double ventaka_omega{};
        /// 分布函数
        DistributionFunction g_t{}, g_bp{};
        DistributionFunction h_t{}, h_bp{};
        std::vector<Vec3D> slope_g{};
        std::vector<Vec3D> slope_h{};
        /// 宏观量
        PhysicalVar::MacroVars macro_vars{};
        /// 构造函数
        explicit Cell(MESH::Cell<int> &cell, Scheme &_solver);
        /// 算法函数
        void get_f_bp();
        void get_grad_f_bp();
        void get_macro_var();
        void update_f_t();
        void update_geom();
    };
    /// Scheme Face
    class Face {
    public:
        Scheme &solver;
        MESH::Face<int> &mesh_face;
        /// 分布函数
        DistributionFunction g{}, g_b{};
        DistributionFunction h{}, h_b{};
        /// 宏观量
        PhysicalVar::MacroVars macro_vars{};
        /// 构造函数
        explicit Face(MESH::Face<int> &face, Scheme &_solver);
        /// 算法函数
        void get_f_b();
        void get_macro_var();
        void get_f();
        void do_boundary();
    };
    /// Container
    std::vector<Cell> CELLS;
    std::vector<Face> FACES;
    inline Cell & get_cell(const int &_key);
    inline Face & get_face(const int &_key);

    /// Solver function
    void init();
    void info() const;
    void do_step();
    void do_save();
    void do_residual();
    void do_crashed(Scheme::Cell &cell);
    void do_crashed(Scheme::Face &face);
};

#endif //SOLVER_DUGKS_SHAKHOV
