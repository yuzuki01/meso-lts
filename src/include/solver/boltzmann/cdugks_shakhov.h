#ifndef SOLVER_CDUGKS_SHAKHOV
#define SOLVER_CDUGKS_SHAKHOV

class CDUGKS_SHAKHOV : public BasicSolver {
public:
    /// Constructor
    using Scheme = CDUGKS_SHAKHOV;
    explicit CDUGKS_SHAKHOV(ConfigReader &_config, ArgParser &_parser);
    /// Check Point
    CheckPoint<CDUGKS_SHAKHOV> check_point;
    /// Mesh
    MESH::StaticMesh phy_mesh{MESH_TYPE_NORMAL, config.phy_mesh};
    MESH::StaticMesh dvs_mesh{MESH_TYPE_NO_FACE, config.dvs_mesh};
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
        MESH::Cell &mesh_cell;
        /// 最小二乘法
        LeastSquare lsp;
        /// Limiter
        double ventaka_omega{};
        /// 分布函数
        DistributionFunction g{};
        DistributionFunction h{};
        std::vector<Vec3D> slope_g{};
        std::vector<Vec3D> slope_h{};
        /// 宏观量
        PhysicalVar::MacroVars macro_vars{};
        /// 构造函数
        explicit Cell(MESH::Cell &cell, Scheme &_solver);
        /// 算法函数
        void init(PhysicalVar::MacroVars init_var);
        void get_grad_f();
        void get_macro_var();
        void update_f();
        void update_geom();
    };
    /// Scheme Face
    class Face {
    public:
        Scheme &solver;
        MESH::Face &mesh_face;
        /// 分布函数
        DistributionFunction g{};
        DistributionFunction h{};
        /// 宏观量
        PhysicalVar::MacroVars macro_vars{};
        /// 构造函数
        explicit Face(MESH::Face &face, Scheme &_solver);
        /// 算法函数
        void get_f_b();
        void get_macro_var();
        void get_f();
        void do_boundary();
    };
    /// Container
    std::vector<Cell> CELLS;
    std::vector<Face> FACES;
    Cell & get_cell(const int &_key);
    Face & get_face(const int &_key);

    /// Solver function
    void init();
    void info() const;
    void do_step();
    void do_save();
    void do_residual();
    void do_crashed(Scheme::Cell &cell);
    void do_crashed(Scheme::Face &face);
};

#endif //SOLVER_CDUGKS_SHAKHOV
