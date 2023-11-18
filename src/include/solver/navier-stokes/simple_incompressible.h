#ifndef SOLVER_SIMPLE_INCOMPRESSIBLE
#define SOLVER_SIMPLE_INCOMPRESSIBLE

class SIMPLE {
public:
    ConfigReader &config;
    ArgParser &parser;
    Logger logger;

    bool is_crashed, continue_to_run = false;
    int step, save_interval, max_step;

    /// Constructor
    using Scheme = SIMPLE;
    explicit SIMPLE(ConfigReader &_config, ArgParser &_parser);
    /// ListMap
    //MESH::ListMesh phy_mesh = MESH::ListMesh();
};

#endif // SOLVER_SIMPLE_INCOMPRESSIBLE
