# meso/mesh

meso 求解器

## 创建自己的求解器

### 1.继承基类

为了统一接口，自定义求解器需要通过继承求解器基类(BasicSolver)，重载基类中的函数来实现。

```c++
class BasicSolver {
public:
    // 日志输出
    Logger logger;
    // 配置文件
    ConfigReader &config;
    // 启动参数
    ArgParser &parser;
    // 通用参数
    int step{}, max_step{}, save_interval{};
    bool is_crashed{}, continue_to_run{}; 
    double stop_at_specific_time{};
    double simulate_time{}, stop_time{};

    explicit BasicSolver(ConfigReader &_config, ArgParser &_parser);
    void init();
    void info();
    void do_save();
    void do_step();
};
```

### 2.头文件声明

创建 MySolver 求解器的头文件。

```c++
/// src/include/solver/boltzmann(PDE to be solved)/MySolver.h

#ifndef SOLVER_MYSOLVER
#define SOLVER_MYSOLVER

class MySolver : public BasicSolver {
public:
    /// Constructor
    using Scheme = MySolver;
    explicit MySolver(ConfigReader &_config, ArgParser &_parser);
    /// Mesh
    MESH::ListMesh mesh = MESH::ListMesh(MESH_TYPE_NORMAL, config.phy_mesh);
    /// Physical
    double Re, Ma;  // any other param you need
    double R, T, Rho, L;
    // Physical Formula
    /// defined here
    
    /// Scheme Cell
    class Cell {
    public:
        Scheme &solver;
        /// other codes
    };
    /// Scheme Face
    class Face {
    public:
        /// other codes
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

#endif // SOLVER_MYSOLVER

```
