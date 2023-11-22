#ifndef SOLVER_GKS
#define SOLVER_GKS

class GKS {
public:
    ConfigReader &config;
    ArgParser &parser;
    Logger logger;

    bool is_crashed, continue_to_run = false;
    int step, save_interval, max_step;

    /// Constructor
    using Scheme = GKS;
    explicit GKS(ConfigReader &_config, ArgParser &_parser);
    /// Mesh
    MESH::ListMesh mesh = MESH::ListMesh(MESH_TYPE_NORMAL, config.phy_mesh);
    /// Physical
    double Ma;
    double R, T, Rho, L;
    double CFL{}, dt{};
    int D;
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

    /// Physical Formula
    void P2C(PhyVar &_var) const; // primitive -> conserved
    void C2P(PhyVar &_var) const; // conserved -> primitive

    /// Scheme Cell
    class Cell {

    };
};


#endif //SOLVER_GKS
