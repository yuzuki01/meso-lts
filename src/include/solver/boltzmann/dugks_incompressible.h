#ifndef SOLVER_DUGKS_INCOMPRESSIBLE
#define SOLVER_DUGKS_INCOMPRESSIBLE

class DUGKS_INCOMPRESSIBLE : public BasicSolver {
public:
    /// Constructor
    using Scheme = DUGKS_INCOMPRESSIBLE;
    explicit DUGKS_INCOMPRESSIBLE(ConfigReader &_config, ArgParser &_parser);
    /// Check Point
    CheckPoint<DUGKS_INCOMPRESSIBLE> check_point;
    /// Mesh
    MESH::StaticMesh phy_mesh{MESH_TYPE_NORMAL, config.phy_mesh};
    MESH::StaticMesh dvs_mesh{MESH_TYPE_NO_FACE, config.dvs_mesh};
    /// Physical
    double Re, Ma;
    double R, T0, Rho0, L;
    double RT{}, CFL{}, dt{}, half_dt{};
    double tau{};
    int D;

    /// Physical Formula
    inline double f_maxwell(double density, const Vec3D &particle_velocity, const Vec3D &flow_velocity) const;
    inline double f_maxwell(const PhysicalVar::MacroVars &macro_var, const Vec3D &particle_velocity) const;

    using DistributionFunction = std::vector<double>;
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell &mesh_cell;
        /// 最小二乘法
        LeastSquare lsp;
        /// 分布函数
        DistributionFunction f_t{}, f_bp{};
        std::vector<Vec3D> slope_f{};
        /// 宏观量
        PhysicalVar::MacroVars macro_vars{};
        /// 构造函数
        explicit Cell(MESH::Cell &cell, Scheme &_solver);
        /// 算法函数
        void init(const PhysicalVar::MacroVars &init_var);
        void get_f_bp();
        void get_grad_f_bp();
        void get_macro_var();
        void update_f_t();
        void update_least_square();
    };
    /// Scheme Face
    class Face {
    public:
        Scheme &solver;
        MESH::Face &mesh_face;
        /// 分布函数
        DistributionFunction f{}, f_b{};
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
};

#endif  // SOLVER_DUGKS_INCOMPRESSIBLE
