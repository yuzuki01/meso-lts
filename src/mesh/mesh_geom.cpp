#include "mesh.h"


/// 局部函数
TP_key
double geom_edge_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_tria_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_quad_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_tria_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_quad_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_tetra_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_wedge_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_pyram_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);

TP_key
double geom_brick_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh);


/// 函数实现
int GEOM::node_num(const int elem_type) {
    switch (elem_type) {
        case EDGE:
            return 2;
        case QUAD:
        case TETRA:
            return 4;
        case TRIA:
            return 3;
        case BRICK:
            return 8;
        case WEDGE:
            return 6;
        case PYRAM:
            return 5;
        default:
            std::string s = "GEOM::node_num() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

int GEOM::face_num(const int elem_type) {
    switch (elem_type) {
        case QUAD:
        case TETRA:
            return 4;
        case TRIA:
            return 3;
        case BRICK:
            return 6;
        case WEDGE:
        case PYRAM:
            return 5;
        default:
            std::string s = "GEOM::face_num() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

double GEOM::vector_angles_2d(const Vec3D &_vec1, const Vec3D &_vec2) {
    double m1 = _vec1.magnitude(), m2 = _vec2.magnitude();
    return acos((_vec1 * _vec2) / (m1 * m2));
}

TP_func
void GEOM::face_normal_vector(MESH::Face<int> &face, MESH::Mesh<int> &mesh) {
    Vec3D vcf = face.position - mesh.get_cell(face.on_cell_key).position;
    if (mesh.dimension() == 2) {
        /// 2D config
        /// get vector along the edge
        Vec3D vfn = mesh.get_node(face.node_key[0]).position - mesh.get_node(face.node_key[1]).position;
        vfn.norm();
        /// get angle
        double angle_x = vector_angles_2d(vfn, {1.0, 0.0, 0.0});
        double angle_y = vector_angles_2d(vfn, {0.0, 1.0, 0.0});
        if (angle_y <= Pi_over_2) {
            angle_x = Pi_over_2 + angle_x;
        } else {
            angle_x = Pi_over_2 - angle_x;
        }
        /// get normal vector
        Vec3D nv(cos(angle_x), sin(angle_x), 0.0);
        /// direction
        if (nv * vcf >= 0.0) {
            face.on_cell_nv = nv;
            face.inv_cell_nv = -nv;
        } else {
            face.on_cell_nv = -nv;
            face.inv_cell_nv = nv;
        }
    } else if (mesh.dimension() == 3) {
        /// 3D config
        /// get near edges on face
        Vec3D v01 = mesh.get_node(face.node_key[1]).position - mesh.get_node(face.node_key[0]).position,
                v12 = mesh.get_node(face.node_key[2]).position - mesh.get_node(face.node_key[1]).position;
        /// get normal vector
        Vec3D nv = v01 ^ v12;   // cross
        nv.norm();
        /// direction
        if (nv * vcf >= 0.0) {
            face.on_cell_nv = nv;
            face.inv_cell_nv = -nv;
        } else {
            face.on_cell_nv = -nv;
            face.inv_cell_nv = nv;
        }
    } else {
        std::string s = "GEOM::face_normal_vector() caught invalid dimension.";
        warn_println(s);
        throw std::invalid_argument(s);
    }
}

TP_func
void GEOM::face_normal_vector(MESH::Face<std::string> &face, MESH::Mesh<std::string> &mesh) {
    Vec3D vcf = face.position - mesh.get_cell(face.on_cell_key).position;
    if (mesh.dimension() == 2) {
        /// 2D config
        /// get vector along the edge
        Vec3D vfn = mesh.get_node(face.node_key[0]).position - mesh.get_node(face.node_key[1]).position;
        vfn.norm();
        /// get angle
        double angle_x = vector_angles_2d(vfn, {1.0, 0.0, 0.0});
        double angle_y = vector_angles_2d(vfn, {0.0, 1.0, 0.0});
        if (angle_y <= Pi_over_2) {
            angle_x = Pi_over_2 + angle_x;
        } else {
            angle_x = Pi_over_2 - angle_x;
        }
        /// get normal vector
        Vec3D nv(cos(angle_x), sin(angle_x), 0.0);
        /// direction
        if (nv * vcf >= 0.0) {
            face.on_cell_nv = nv;
            face.inv_cell_nv = -nv;
        } else {
            face.on_cell_nv = -nv;
            face.inv_cell_nv = nv;
        }
    } else if (mesh.dimension() == 3) {
        /// 3D config
        /// get near edges on face
        Vec3D v01 = mesh.get_node(face.node_key[1]).position - mesh.get_node(face.node_key[0]).position,
                v12 = mesh.get_node(face.node_key[2]).position - mesh.get_node(face.node_key[1]).position;
        /// get normal vector
        Vec3D nv = v01 ^ v12;   // cross
        nv.norm();
        /// direction
        if (nv * vcf >= 0.0) {
            face.on_cell_nv = nv;
            face.inv_cell_nv = -nv;
        } else {
            face.on_cell_nv = -nv;
            face.inv_cell_nv = nv;
        }
    } else {
        std::string s = "GEOM::face_normal_vector() caught invalid dimension.";
        warn_println(s);
        throw std::invalid_argument(s);
    }
}

TP_func
Vec3D GEOM::face_position(const MESH::Face<int> &face, MESH::Mesh<int> &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : face.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(face.node_key.size());
}

TP_func
Vec3D GEOM::face_position(const MESH::Face<std::string> &face, MESH::Mesh<std::string> &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : face.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(face.node_key.size());
}

TP_func
Vec3D GEOM::cell_position(const MESH::Cell<int> &cell, MESH::Mesh<int> &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : cell.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(cell.node_key.size());
}

TP_func
Vec3D GEOM::cell_position(const MESH::Cell<std::string> &cell, MESH::Mesh<std::string> &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : cell.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(cell.node_key.size());
}

TP_func
double GEOM::face_area(const MESH::Face<int> &face, MESH::Mesh<int> &mesh) {
    switch (face.type) {
        case EDGE:
            return geom_edge_area(face.node_key, mesh);
        case TETRA:
            return geom_tria_area(face.node_key, mesh);
        case QUAD:
            return geom_quad_area(face.node_key, mesh);
        default:
            std::string s = "GEOM::face_area() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::face_area(const MESH::Face<std::string> &face, MESH::Mesh<std::string> &mesh) {
    switch (face.type) {
        case EDGE:
            return geom_edge_area(face.node_key, mesh);
        case TETRA:
            return geom_tria_area(face.node_key, mesh);
        case QUAD:
            return geom_quad_area(face.node_key, mesh);
        default:
            std::string s = "GEOM::face_area() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::cell_volume(const MESH::Cell<int> &cell, MESH::Mesh<int> &mesh) {
    switch (cell.type) {
        case TRIA:
            return geom_tria_volume(cell.node_key, mesh);
        case QUAD:
            return geom_quad_volume(cell.node_key, mesh);
        case TETRA:
            return geom_tetra_volume(cell.node_key, mesh);
        case PYRAM:
            return geom_pyram_volume(cell.node_key, mesh);
        case WEDGE:
            return geom_wedge_volume(cell.node_key, mesh);
        case BRICK:
            return geom_brick_volume(cell.node_key, mesh);
        default:
            std::string s = "GEOM::cell_volume() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::cell_volume(const MESH::Cell<std::string> &cell, MESH::Mesh<std::string> &mesh) {
    switch (cell.type) {
        case TRIA:
            return geom_tria_volume(cell.node_key, mesh);
        case QUAD:
            return geom_quad_volume(cell.node_key, mesh);
        case TETRA:
            return geom_tetra_volume(cell.node_key, mesh);
        case PYRAM:
            return geom_pyram_volume(cell.node_key, mesh);
        case WEDGE:
            return geom_wedge_volume(cell.node_key, mesh);
        case BRICK:
            return geom_brick_volume(cell.node_key, mesh);
        default:
            std::string s = "GEOM::cell_volume() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

/// 局部函数实现
TP_key
double geom_edge_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]);
    return (n0.position - n1.position).magnitude();
}

TP_key
double geom_tria_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]);
    Vec3D v01 = n1.position - n0.position, v02 = n2.position - n0.position;
    return (v01 ^ v02).magnitude() / 2.0;
}

TP_key
double geom_quad_area(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    return geom_tria_area({node_keys[0], node_keys[1], node_keys[2]}, mesh)
           + geom_tria_area({node_keys[0], node_keys[3], node_keys[2]}, mesh);
}

TP_key
double geom_tria_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]);
    Vec3D v01 = n1.position - n0.position, v02 = n2.position - n0.position;
    return (v01 ^ v02).magnitude() / 2.0;
}

TP_key
double geom_quad_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    return geom_tria_volume({node_keys[0], node_keys[1], node_keys[3]}, mesh)
           + geom_tria_volume({node_keys[2], node_keys[1], node_keys[3]}, mesh);
}

TP_key
double geom_tetra_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    /// top point - 0
    /// bottom tria - 1, 2, 3
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]),
            &n3 = mesh.get_node(node_keys[3]);
    Vec3D v01 = n1.position - n0.position,
            v02 = n2.position - n0.position,
            v03 = mesh.get_node(node_keys[3]).position - n0.position;
    return fabs(v01 * (v02 ^ v03)) / 6.0;
}

TP_key
double geom_pyram_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// basic shape
    /// top point - 4
    /// bottom quad - 0, 1, 3, 2
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]),
            &n3 = mesh.get_node(node_keys[3]),
            &n4 = mesh.get_node(node_keys[4]);
    Vec3D v01 = n1.position - n0.position,
            v02 = n2.position - n0.position,
            v04 = n4.position - n0.position,
            v31 = n1.position - n3.position,
            v32 = n2.position - n3.position;
    return fabs(v04 * (v01 ^ v02)) / 6. + fabs(v04 * (v31 ^ v32)) / 6.;
}

TP_key
double geom_wedge_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// top - tria 0, 1, 2
    /// bottom - tria 3, 4, 5
    return geom_pyram_volume({node_keys[0], node_keys[2], node_keys[3], node_keys[5], node_keys[1]}, mesh)
           + geom_tetra_volume({node_keys[1], node_keys[3], node_keys[4], node_keys[5]}, mesh);
}

TP_key
double geom_brick_volume(const std::vector<key_type> &node_keys, MESH::Mesh<key_type> &mesh) {
    /// wedge 02-13-46
    double wedge014236 = geom_wedge_volume({
                                                   node_keys[0], node_keys[1], node_keys[4],
                                                   node_keys[2], node_keys[3], node_keys[6]
                                           }, mesh);
    /// wedge 57-13-46
    double wedge145367 = geom_wedge_volume({
                                                   node_keys[1], node_keys[4], node_keys[5],
                                                   node_keys[3], node_keys[6], node_keys[7]
                                           }, mesh);
    return wedge014236 + wedge145367;
}
