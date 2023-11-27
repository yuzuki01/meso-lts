#ifndef SOLVER_GKS
#define SOLVER_GKS

class GKS : public BasicSolver {
private:
    const std::unordered_map<std::string, int> gks_solver_map{
            {"kfvs1st", 0}, {"kfvs2nd", 1}
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
    double R, T, Rho, L;
    double CFL;
    int D, K, stage, gks_solver_key;
    double c1_euler{}, c2_euler{}, dt{}, gamma{};
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
    struct PhyVar_lambda {
        double density = 0.0;                        // prime[0]
        Vec3D velocity{0.0, 0.0, 0.0};     // prime[1-3]
        double lambda = 0.0;                         // prime[4]
    };
    struct AVector {
        double a1;
        Vec3D a234;
        double a5;
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
    double time_coefficient[5][5][3]{};

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
        double upvxi[7][7][7][3];
        double unvxi[7][7][7][3];
        double uvxi[7][7][7][3];
        double xi2;
        double xi4;
        MMDF();
        MMDF(const Vec3D &velocity_in, double lambda_in);
        void calcualte_MMDF();
    };
    inline void Prim_to_Cons(PhyVar &_var) const; // primitive -> conserved
    inline void Cons_to_Prim(PhyVar &_var) const; // conserved -> primitive
    inline void Cons_to_Char(PhyVar &_var) const;
    inline void Char_to_Cons(PhyVar &_var) const;
    inline void Convar_to_ULambda(double* _prim, PhyVar &_var) const; //change conservative variables to rho u lambda

    ///GKS calculate tau
    double Get_Tau_NS(double density0, double lambda0) const;
    double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt) const;//cacluate tau
    double TauNS_Sutherland(double density0, double lambda0);
    double TauNS_power_law(double density0, double lambda0);
    /// GKS collision
    void Collision(PhyVar &center, double _left, double _right, MMDF& m2, MMDF& m3);
    ///solution of matrix equation b=Ma 3D
    void A_point(double* a, const std::vector<double> &der1i, double* prim);
    //a general G function
    /**
     * (4,11)
     * 1 0 0 0 for Wc
     * 0 0 0 0 for partialW
     **/
    void G_address(int no_u, int no_v,int no_w, int no_xi, double* psi, double a[5], MMDF& m);
    void GL_address(int no_u, int no_v, int no_w, int no_xi, double* psi, double a[4], MMDF& m);
    void GR_address(int no_u, int no_v, int no_w, int no_xi, double* psi, double a[4], MMDF& m);
    /// GKS solvers
    const int kfvs1st = gks_solver_map.at("kfvs1st");

    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        MESH::Cell<int> *mesh_cell_ptr;
        /// Geom
        Mat3D coordinate_trans;
        LeastSquare lsp;

        PhyVar pv;

        explicit Cell(MESH::Cell<int> &cell, Scheme &_solver);
        void update(int now_stage);
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

        void reconstruct(int now_stage);
        void calculate_flux(int now_stage);
        void do_boundary(int now_stage);
    };
    /// Container
    std::vector<Cell> CELLS, GHOST_CELLS;
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
