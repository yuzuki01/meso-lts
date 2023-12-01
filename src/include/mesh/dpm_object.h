namespace DPM {
    TP_mesh
    class Particle;
}

TP_mesh
class DPM::Particle {
public:
    int key{};
    double diameter{};
    Vec3D position{0.0, 0.0, 0.0};
    mesh_type &mesh;     // const
    int in_cell_key{};
    MESH::Cell *in_cell_ptr{};

    Particle<mesh_type>(int _key, double _diameter, const Vec3D &_position,
                       mesh_type &_mesh, MESH::Cell &_in_cell) : key(_key), diameter(_diameter),
                                                                           position(_position),
                                                                           mesh(_mesh),
                                                                           in_cell_key(_in_cell.key),
                                                                           in_cell_ptr(&_in_cell) {};

    bool is_in_cell(int cell_key);

    /// code in dpm_geom.cpp
    bool is_in_cell(MESH::Cell &cell);
};
