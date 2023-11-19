#ifndef SOLVER_SIMPLE_INCOMPRESSIBLE
#define SOLVER_SIMPLE_INCOMPRESSIBLE

class SIMPLE_INCOMPRESSIBLE {
public:
    ConfigReader &config;
    ArgParser &parser;
    Logger logger;

    bool is_crashed, continue_to_run = false;
    int step, save_interval, max_step;

    /// Constructor
    using Scheme = SIMPLE_INCOMPRESSIBLE;
    explicit SIMPLE_INCOMPRESSIBLE(ConfigReader &_config, ArgParser &_parser);
    /// ListMap
    MESH::ListMesh mesh = MESH::ListMesh(MESH_TYPE_NORMAL, "mesh");
    /// Physical
    double Re, Ma;
    double Rho, Length;
    double CFL{}, dt{};
    /// Physical Vars
    struct PhysicalVars {
        double density{};
        double pressure{};
        Vec3D velocity{};
    };
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;
        /// 最小二乘法
        LeastSquare lsp;
        /// 物理量
        PhysicalVars phy_vars{};
        /// 构造函数
        explicit Cell(MESH::ListMesh &cell, Scheme &_solver);
        /// 算法函数

    };
    /// Scheme Face
    class Face {
    public:
        Scheme &solver;
        MESH::Face<int> *mesh_face_ptr;
        /// 物理量
        PhysicalVars phy_var{};
        /// 构造函数
        explicit Face(MESH::ListMesh &cell, Scheme &_solver);
        /// 算法函数

    };
    /// Container
    std::vector<Cell> CELLS;
    std::vector<Face> FACES;
    Cell & get_cell(const int &_key);
    Face & get_face(const int &_key);

    /// Solver function
    void init();
    void info();
    void do_step();
    void do_save();
    void do_residual();
    void do_crashed(Scheme::Cell &cell);
};

#endif // SOLVER_SIMPLE_INCOMPRESSIBLE
