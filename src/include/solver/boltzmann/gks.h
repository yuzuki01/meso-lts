#ifndef SOLVER_GKS
#define SOLVER_GKS

class GKS : public BasicSolver {
public:
    using Scheme = GKS;
    explicit GKS(ConfigReader &_config, ArgParser &_parser);
    CheckPoint<Scheme> check_point;
    MESH::StaticMesh mesh{MeshTypeNormal, config.phy_mesh};
    double Ma;
    double R, T0, Rho0, L, gamma;
    double CFL, dt{};
    int D, K;

    /// Physical Formula
    inline double lambda(double temperature) const;
    class MMDF {
    private:
        int zone, K;
        bool is_computed[7]{};
        double _moment[7]{};
    public:
        enum {full_zone, left_zone, right_zone};
        double u, lambda;
        /// moment of maxwell distribution function
        MMDF(double _u, double _lambda, int _k,  int integral_zone);
        double operator()(int _o);          // macro
        double operator[](int _o) const;    // internal
    };

    class Cell {
    public:
        Scheme &solver;
        MESH::Cell &mesh_cell;
        LeastSquare lsp;
        Physical::MacroVars macro_vars{};
        explicit Cell(MESH::Cell &cell, Scheme &_solver);
        void init(const Physical::MacroVars &init_var);
        void update();
    };

    class Face {
    public:
        Scheme &solver;
        MESH::Face &mesh_face;
        Physical::MacroVars w_left, w_right, w_c;
        Physical::ConservedFlux flux;
        explicit Face(MESH::Face &face, Scheme &_solver);
        void boundary();
        void reconstruct();
    };

    std::vector<Cell> CELLS;
    std::vector<Face> FACES;
    Cell & get_cell(const int &_key);
    Face & get_face(const int &_key);

    void init();
    void info() const;
    void do_step();
    void do_save();
    void do_residual();
};


/// check_point

TP_func void CheckPoint<GKS>::init_field(const Physical::MacroVars &_var);

TP_func void CheckPoint<GKS>::init_from_file(const std::string &file_path);

TP_func void CheckPoint<GKS>::write_to_file(const std::string &file_path);


#endif //SOLVER_GKS
