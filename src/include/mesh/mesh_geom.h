/**
 * included by mesh.h
 */

namespace GEOM {
    /**
     * 1 = Edge             线段
     * 2 = Quadrilateral    四边形
     * 3 = Triangle         三角形
     * 4 = Brick            立方体
     * 5 = Wedge (Prism)    三棱柱
     * 6 = Tetrahedron      三棱锥
     * 7 = Pyramid          四棱锥
     */
    const int EDGE = 1;
    const int QUAD = 2;
    const int TRIA = 3;
    const int BRICK = 4;
    const int WEDGE = 5;
    const int TETRA = 6;
    const int PYRAM = 7;

    /// get number of node/face
    int node_num(int elem_type);
    int face_num(int elem_type);
    /// calculation
    TP_key
    void face_normal_vector(MESH::Face<key_type> &face, MESH::Mesh<key_type> &mesh);
    TP_key
    Vec3D face_position(const MESH::Face<key_type> &face, MESH::Mesh<key_type> &mesh);
    TP_key
    Vec3D cell_position(const MESH::Cell<key_type> &cell, MESH::Mesh<key_type> &mesh);
    TP_key
    double face_area(const MESH::Face<key_type> &face, MESH::Mesh<key_type> &mesh);
    TP_key
    double cell_volume(const MESH::Cell<key_type> &cell, MESH::Mesh<key_type> &mesh);
    double vector_angles_2d(const Vec3D &_vec1, const Vec3D &_vec2);
    /// particle tracker
    TP_key
    bool is_particle_in_cell(MESH::Cell<key_type> &cell, MESH::Mesh<key_type> &mesh);
}
