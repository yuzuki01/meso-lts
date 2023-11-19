namespace DPM {
    TP_key_mesh
    class Particle;
}

TP_key_mesh
class DPM::Particle {
public:
    key_type key{};
    double diameter{};
    Vec3D position{0.0, 0.0, 0.0};
    mesh_type &mesh;     // const
    key_type in_cell_key{};
    MESH::Cell<key_type> *in_cell_ptr{};

    Particle<key_type, mesh_type>(key_type _key, double _diameter, const Vec3D &_position,
                       mesh_type &_mesh, MESH::Cell<key_type> &_in_cell) : key(_key), diameter(_diameter),
                                                                           position(_position),
                                                                           mesh(_mesh),
                                                                           in_cell_key(_in_cell.key),
                                                                           in_cell_ptr(&_in_cell) {};

    bool is_in_cell(const key_type &cell_key);

    /// code in dpm_geom.cpp
    bool is_in_cell(MESH::Cell<key_type> &cell);
};
