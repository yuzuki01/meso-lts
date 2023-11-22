#ifndef SOLVER_GKS
#define SOLVER_GKS

class GKS : public BasicSolver {
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
    int D, K, stage;
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
    };
    struct GradVar {
        Vec3D density{0.0, 0.0, 0.0};
        Vec3D energy{0.0, 0.0, 0.0};
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
    void P2C(PhyVar &_var) const; // primitive -> conserved
    void C2P(PhyVar &_var) const; // conserved -> primitive

    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;

        int cell_stage;
        PhyVar pv;
        void update();
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
