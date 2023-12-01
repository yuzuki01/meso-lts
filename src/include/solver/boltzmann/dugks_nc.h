#ifndef SOLVER_DUGKS_NC
#define SOLVER_DUGKS_NC

class DUGKS_NC : public BasicSolver {
public:
    /// Constructor
    using Scheme = DUGKS_NC;
    explicit DUGKS_NC(ConfigReader &_config, ArgParser &_parser);
    /// Mesh
    MESH::StaticMesh phy_mesh{MESH_TYPE_NORMAL, config.phy_mesh};
    MESH::StaticMesh dvs_mesh{MESH_TYPE_NO_FACE, "GaussHermit"};
    /// Physical
    double Re, Ma, Ra, Pr;
    double R, T0, dT, Rho0, L;
    double RT{}, CFL{}, dt{}, half_dt{};
    double tau_f{}, tau_g{}, beta_g{};
    int D;
    /// Macro Physical Variables
    struct MacroVars {
        double density = 0.0;
        double temperature = 0.0;
        Vec3D velocity{0.0, 0.0, 0.0};
        Vec3D heat_flux{0.0, 0.0, 0.0};
    };

    /// Physical Formula
    inline double f_maxwell(double density, const Vec3D &particle_velocity, const Vec3D &flow_velocity) const;
    inline double g_maxwell(double temperature, const Vec3D &particle_velocity, const Vec3D &);

    using DistributionFunction = std::vector<double>;
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;
        /// 最小二乘法
        LeastSquare lsp;
        /// 分布函数
        DistributionFunction f_t{}, f_bp{};
        std::vector<Vec3D> slope_f{};
        /// 宏观量
        MacroVars macro_vars{};
        /// 构造函数
        explicit Cell(MESH::Cell<int> &cell, Scheme &_solver);
        /// 算法函数
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
        MESH::Face<int> *mesh_face_ptr;
        /// 分布函数
        DistributionFunction f{}, f_b{};
        /// 宏观量
        MacroVars macro_vars{};
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
    Cell & get_cell(const int &_key);
    Face & get_face(const int &_key);

    /// Solver function
    void init();
    void info() const;
    void do_step();
    void do_save();
};

#endif  // SOLVER_DUGKS_NC
