/**
 * included by solver.h
 */

class ConfigReader : BasicReader {
protected:
    /// config file mark
    const std::string def_mark = "def:";
    const std::string bc_mark = "boundary:";
    const std::string end_mark = ":end";
public:
    using ParamMap = std::unordered_map<std::string, std::string>;
    std::string name{}, phy_mesh{}, dvs_mesh{}, solver{};
    int thread = 1;
    ParamMap param{};
    std::vector<BoundaryParam> mark_sets{};
    ConfigReader() = default;
    explicit ConfigReader(const std::string &file_path);
    void info();

    std::string operator[](const std::string &key);
    template<class T> T get(const std::string &_key);

    bool set_mesh_mark(MESH::StaticMesh &mesh);
    void set_mesh_mark(MESH::MapMesh &mesh);
};
