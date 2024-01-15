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
    MESH::StaticMesh phy_mesh{MeshTypeNormal, config.phy_mesh};
    MESH::StaticMesh dvs_mesh{MeshTypeNoFace, config.dvs_mesh};
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
        DistributionFunction g_t{}, g_bp{};
        DistributionFunction h_t{}, h_bp{};
        std::vector<Vec3D> slope_g{};
        std::vector<Vec3D> slope_h{};
        /// 宏观量
        Physical::MacroVars macro_vars{};
        /// 构造函数
        explicit Cell(MESH::Cell &cell, Scheme &_solver);
        /// 算法函数
        void init(Physical::MacroVars init_var);
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
        MESH::Face &mesh_face;
        /// 分布函数
        DistributionFunction g{}, g_b{};
        DistributionFunction h{}, h_b{};
        /// 宏观量
        Physical::MacroVars macro_vars{};
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

/// check point

TP_func void CheckPoint<DUGKS_SHAKHOV>::init_field(const Physical::MacroVars &_var);

TP_func void CheckPoint<DUGKS_SHAKHOV>::init_from_file(const std::string &file_path);

TP_func void CheckPoint<DUGKS_SHAKHOV>::write_to_file(const std::string &file_path);

#endif //SOLVER_DUGKS_SHAKHOV
