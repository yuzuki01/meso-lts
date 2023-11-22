#ifndef SOLVER_DUGKS_SHAKHOV
#define SOLVER_DUGKS_SHAKHOV

class DUGKS_SHAKHOV : public BasicSolver {
public:
    /// Constructor
    using Scheme = DUGKS_SHAKHOV;
    explicit DUGKS_SHAKHOV(ConfigReader &_config, ArgParser &_parser);
    /// Mesh
    MESH::ListMesh phy_mesh{MESH_TYPE_NORMAL, config.phy_mesh};
    MESH::ListMesh dvs_mesh{MESH_TYPE_NO_FACE, config.dvs_mesh};
    /// Physical
    double Re, Ma;
    double R, T, Rho, L;
    double RT{}, CFL{}, dt{}, half_dt{};
    double tau{};
    int D, K;
    /// Macro Physical Variables
    struct MacroVars {
        double density = 0.0;
        double temperature = 0.0;
        double pressure = 0.0;
        Vec3D velocity{0.0, 0.0, 0.0};
        Vec3D heat_flux{0.0, 0.0, 0.0};
    };

    /// Physical Formula
    inline double g_maxwell(double density, double temperature,
                            const Vec3D &macro_velocity, const Vec3D &particle_velocity) const;
    inline double g_shakhov(double g_m, const Vec3D &particle_velocity, const Vec3D heat_flux) const;
    inline double h_maxwell(double g_m, double temperature) const;
    inline double h_shakhov(double g_m, const Vec3D &particle_velocity, const Vec3D heat_flux) const;

    using DistributionFunction = std::vector<double>;
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;
        /// 最小二乘法
        LeastSquare lsp;
        /// 分布函数
        DistributionFunction g_t{}, g_bp{};
        DistributionFunction h_t{}, h_bp{};
        std::vector<Vec3D> slope_g{};
        std::vector<Vec3D> slope_h{};
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
        DistributionFunction g{}, g_b{};
        DistributionFunction h{}, h_b{};
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

#endif //SOLVER_DUGKS_SHAKHOV
