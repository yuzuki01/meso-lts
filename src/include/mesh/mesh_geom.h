/**
 * included by mesh.h
 */

namespace GEOM {
    /**
     * su2 geom id
     **/
    const int LINE = 3;         // 线
    const int TRIA = 5;         // 三角形
    const int QUAD = 9;         // 四边形
    const int TETRA = 10;       // 四面体
    const int HEXAH = 12;       // 六面体
    const int PRISM = 13;       // 金字塔
    const int PYRAM = 14;       // 三棱柱

    /// generate face-key for finding in a map
    TP_key
    std::string generate_face_key(key_vector &node_set);
    /// get number of node/face
    int node_num(int elem_type);
    int face_num(int elem_type);
    /// calculation
    TP_key_mesh
    void face_normal_vector(MESH::Face<key_type> &face, mesh_type &mesh);
    TP_key_mesh
    Vec3D face_position(const MESH::Face<key_type> &face, mesh_type &mesh);
    TP_key_mesh
    Vec3D cell_position(const MESH::Cell<key_type> &cell, mesh_type &mesh);
    TP_key_mesh
    double face_area(const MESH::Face<key_type> &face, mesh_type &mesh);
    TP_key_mesh
    double cell_volume(const MESH::Cell<key_type> &cell, mesh_type &mesh);

    double vector_angles_2d(const Vec3D &_vec1, const Vec3D &_vec2);
    /// particle tracker
    TP_key_mesh
    bool is_particle_in_cell(MESH::Cell<key_type> &cell, mesh_type &mesh);
}
