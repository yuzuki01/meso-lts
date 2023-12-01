TP_mesh
class MeshWriter {
private:
    bool is_file_open;
    std::ofstream fp;
    mesh_type &mesh;
    std::string path;
public:
    MeshWriter(std::string _path, mesh_type &_mesh);
    bool is_open() const;
    void write_head(const string_vector &values);
    void write_node();
    void write_data(const std::vector<double> &data);
    void write_geom();
    void close();
};
