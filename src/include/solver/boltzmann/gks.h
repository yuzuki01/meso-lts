#ifndef SOLVER_GKS
#define SOLVER_GKS

class GKS : public BasicSolver {
private:
    const std::unordered_map<std::string, int> gks_solver_map{
            {"kfvs1st", 0},
    };
public:
    /// Constructor
    using Scheme = GKS;
    explicit GKS(ConfigReader &_config, ArgParser &_parser);
    /// Mesh
    MESH::ListMesh mesh = MESH::ListMesh(MESH_TYPE_NORMAL, config.phy_mesh);
    /// Physical
    bool is_prandtle_fix;
    double Ma, Re, Pr;
    double R, T, Rho, L, gamma;
    double c1_euler, c2_euler;
    int D, K, stage, gks_solver_key;
    double CFL{}, dt{};
    /// Physical Variables
    struct PhyVar {
        /// Conserved
        double density = 0.0;
        double energy = 0.0;
        Vec3D momentum{0.0, 0.0, 0.0};
        /// Primitive
        double pressure = 0.0;
        double temperature = 0.0;
        Vec3D velocity{0.0, 0.0, 0.0};
        /// Characteristic
        double char_density = 0.0;
        double char_energy = 0.0;
        Vec3D char_momentum{0.0, 0.0, 0.0};
    };
    struct GradVar {
        Vec3D density{0.0, 0.0, 0.0};
        Vec3D energy{0.0, 0.0, 0.0};
        /// momentum x,y,z
        Vec3D momentum_x{0.0, 0.0, 0.0};
        Vec3D momentum_y{0.0, 0.0, 0.0};
        Vec3D momentum_z{0.0, 0.0, 0.0};
    };
    struct FaceFlux {
        PhyVar flux;
        GradVar der_flux;
    };
    double time_coefficient[5][5][3];

    /// Physical Formula
    class MMDF {
    private:
        Vec3D velocity;
        double lambda;

    public:
        double uwhole[7];
        double uplus[7];
        double uminus[7];
        double vwhole[7];
        double upvxi[7][7][3];
        double unvxi[7][7][3];
        double uvxi[7][7][3];
        double xi2;
        double xi4;
        MMDF();
        MMDF(double u_in, double v_in, double lambda_in);
        void calcualte_MMDF();
    };
    void Prim_to_Cons(PhyVar &_var) const; // primitive -> conserved
    void Cons_to_Prim(PhyVar &_var) const; // conserved -> primitive
    void Cons_to_Char(PhyVar &_var) const;
    void Char_to_Cons(PhyVar &_var) const;

    /// GKS solvers
    const int kfvs1st = gks_solver_map.at("kfvs1st");

    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;
        /// Geom
        double coordinate_trans[3][3];
        LeastSquare lsp;

        int cell_stage;
        PhyVar pv;

        explicit Cell(MESH::Cell<int> &cell, Scheme &_solver);
        void update();
        void set_geom();
    };

    /// Scheme Face
    class Face {
    public:
        Scheme &solver;
        MESH::Face<int> *mesh_face_ptr;
        /// GKS face params
        bool is_reduce_order;
        PhyVar pv_left, pv_right, pv_center;
        GradVar grad_left, grad_right, grad_center;
        FaceFlux fluxes;

        explicit Face(MESH::Face<int> &face, Scheme &_solver);

        void reconstruct();
        void calculate_flux();
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


#endif //SOLVER_GKS
