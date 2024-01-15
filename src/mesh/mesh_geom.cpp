#include "mesh.h"


/// 局部函数
TP_mesh
double geom_line_area(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_tria_area(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_quad_area(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_tria_volume(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_quad_volume(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_tetra_volume(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_prism_volume(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_pyram_volume(const key_vector &node_keys, mesh_type &mesh);

TP_mesh
double geom_hexah_volume(const key_vector &node_keys, mesh_type &mesh);


/// 函数实现
int GEOM::node_num(const int elem_type) {
    switch (elem_type) {
        case LINE:
            return 2;
        case TRIA:
            return 3;
        case QUAD:
        case TETRA:
            return 4;
        case PRISM:
            return 5;
        case PYRAM:
            return 6;
        case HEXAH:
            return 8;
        default:
            std::string s = "GEOM::node_num() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

int GEOM::face_num(const int elem_type) {
    switch (elem_type) {
        case TRIA:
            return 3;
        case QUAD:
        case TETRA:
            return 4;
        case PRISM:
        case PYRAM:
            return 5;
        case HEXAH:
            return 6;
        default:
            std::string s = "GEOM::face_num() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

/// generate key
std::string GEOM::generate_face_key(key_vector &node_set) {
    /// transfer 'string' to 'int' weight.
    const int len = int(node_set.size());
    /// find position where is min.
    int pos = 0;
    auto min = node_set[pos];
    for (int i = 0; i < len; i++) {
        if (node_set[i] < min) {
            pos = i;
            min = node_set[i];
        }
    }
    /// direction
    bool direction;
    int left_pos, right_pos;
    left_pos = (pos == 0) ? (len - 1) : (pos - 1);
    right_pos = (pos == len - 1) ? 0 : (pos + 1);
    direction = node_set[right_pos] < node_set[left_pos];
    /// append key
    int count = 0;
    std::stringstream ss;
    while (count < len) {
        ss << node_set[pos] << ">";
        pos = direction ? (pos + 1) : (pos - 1);
        /// make valid range
        if (pos < 0) pos = len - 1;
        if (pos > len - 1) pos = 0;
        count++;
    }
    return ss.str();
}


double GEOM::vector_angles_2d(const Vec3D &_vec1, const Vec3D &_vec2) {
    double m1 = _vec1.magnitude(), m2 = _vec2.magnitude();
    return acos((_vec1 * _vec2) / (m1 * m2));
}

TP_func
void GEOM::face_normal_vector(MESH::Face &face, MESH::StaticMesh &mesh) {
    Vec3D vcf = face.position - mesh.get_cell(face.on_cell_key).position;
    if (mesh.dimension() == 2) {
        /// 2D config
        /// get vector along the LINE
        Vec3D vfn = mesh.get_node(face.node_key[0]).position - mesh.get_node(face.node_key[1]).position;
        vfn = vfn.norm();
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
        /// get near LINEs on face
        Vec3D v01 = mesh.get_node(face.node_key[1]).position - mesh.get_node(face.node_key[0]).position,
                v12 = mesh.get_node(face.node_key[2]).position - mesh.get_node(face.node_key[1]).position;
        /// get normal vector
        Vec3D nv = v01 ^v12;   // cross
        nv = nv.norm();
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
void GEOM::face_normal_vector(MESH::Face &face, MESH::MapMesh &mesh) {
    Vec3D vcf = face.position - mesh.get_cell(face.on_cell_key).position;
    if (mesh.dimension() == 2) {
        /// 2D config
        /// get vector along the LINE
        Vec3D vfn = mesh.get_node(face.node_key[0]).position - mesh.get_node(face.node_key[1]).position;
        vfn = vfn.norm();
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
        /// get near LINEs on face
        Vec3D v01 = mesh.get_node(face.node_key[1]).position - mesh.get_node(face.node_key[0]).position,
                v12 = mesh.get_node(face.node_key[2]).position - mesh.get_node(face.node_key[1]).position;
        /// get normal vector
        Vec3D nv = v01 ^v12;   // cross
        nv = nv.norm();
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
Vec3D GEOM::face_position(const MESH::Face &face, MESH::StaticMesh &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : face.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(node_num(face.type));
}

TP_func
Vec3D GEOM::face_position(const MESH::Face &face, MESH::MapMesh &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : face.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(node_num(face.type));
}

TP_func
Vec3D GEOM::cell_position(const MESH::Cell &cell, MESH::StaticMesh &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : cell.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(node_num(cell.type));
}

TP_func
Vec3D GEOM::cell_position(const MESH::Cell &cell, MESH::MapMesh &mesh) {
    Vec3D position(0.0, 0.0, 0.0);
    for (auto &node_key : cell.node_key) {
        auto &node = mesh.get_node(node_key);
        position += node.position;
    }
    return position / double(node_num(cell.type));
}

TP_func
double GEOM::face_area(const MESH::Face &face, MESH::StaticMesh &mesh) {
    switch (face.type) {
        case LINE:
            return geom_line_area<MESH::StaticMesh>(face.node_key, mesh);
        case TETRA:
            return geom_tria_area<MESH::StaticMesh>(face.node_key, mesh);
        case QUAD:
            return geom_quad_area<MESH::StaticMesh>(face.node_key, mesh);
        default:
            std::string s = "GEOM::face_area() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::face_area(const MESH::Face &face, MESH::MapMesh &mesh) {
    switch (face.type) {
        case LINE:
            return geom_line_area<MESH::MapMesh>(face.node_key, mesh);
        case TETRA:
            return geom_tria_area<MESH::MapMesh>(face.node_key, mesh);
        case QUAD:
            return geom_quad_area<MESH::MapMesh>(face.node_key, mesh);
        default:
            std::string s = "GEOM::face_area() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::cell_volume(const MESH::Cell &cell, MESH::StaticMesh &mesh) {
    switch (cell.type) {
        case TRIA:
            return geom_tria_volume<MESH::StaticMesh>(cell.node_key, mesh);
        case QUAD:
            return geom_quad_volume<MESH::StaticMesh>(cell.node_key, mesh);
        case TETRA:
            return geom_tetra_volume<MESH::StaticMesh>(cell.node_key, mesh);
        case PYRAM:
            return geom_pyram_volume<MESH::StaticMesh>(cell.node_key, mesh);
        case PRISM:
            return geom_prism_volume<MESH::StaticMesh>(cell.node_key, mesh);
        case HEXAH:
            return geom_hexah_volume<MESH::StaticMesh>(cell.node_key, mesh);
        default:
            std::string s = "GEOM::cell_volume() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

TP_func
double GEOM::cell_volume(const MESH::Cell &cell, MESH::MapMesh &mesh) {
    switch (cell.type) {
        case TRIA:
            return geom_tria_volume<MESH::MapMesh>(cell.node_key, mesh);
        case QUAD:
            return geom_quad_volume<MESH::MapMesh>(cell.node_key, mesh);
        case TETRA:
            return geom_tetra_volume<MESH::MapMesh>(cell.node_key, mesh);
        case PYRAM:
            return geom_pyram_volume<MESH::MapMesh>(cell.node_key, mesh);
        case PRISM:
            return geom_prism_volume<MESH::MapMesh>(cell.node_key, mesh);
        case HEXAH:
            return geom_hexah_volume<MESH::MapMesh>(cell.node_key, mesh);
        default:
            std::string s = "GEOM::cell_volume() caught unsupported geom type.";
            warn_println(s);
            throw std::invalid_argument(s);
    }
}

/// 局部函数实现
TP_mesh
double geom_line_area(const key_vector &node_keys, mesh_type &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]);
    return (n0.position - n1.position).magnitude();
}

TP_mesh
double geom_tria_area(const key_vector &node_keys, mesh_type &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]);
    Vec3D v01 = n1.position - n0.position, v02 = n2.position - n0.position;
    return (v01 ^ v02).magnitude() / 2.0;
}

TP_mesh
double geom_quad_area(const key_vector &node_keys, mesh_type &mesh) {
    return geom_tria_area<mesh_type>({node_keys[0], node_keys[1], node_keys[2]}, mesh)
           + geom_tria_area<mesh_type>({node_keys[0], node_keys[3], node_keys[2]}, mesh);
}

TP_mesh
double geom_tria_volume(const key_vector &node_keys, mesh_type &mesh) {
    /// basic shape
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]);
    Vec3D v01 = n1.position - n0.position, v02 = n2.position - n0.position;
    return (v01 ^ v02).magnitude() / 2.0;
}

TP_mesh
double geom_quad_volume(const key_vector &node_keys, mesh_type &mesh) {
    /// basic shape
    return geom_tria_volume<mesh_type>({node_keys[0], node_keys[1], node_keys[3]}, mesh)
           + geom_tria_volume<mesh_type>({node_keys[2], node_keys[1], node_keys[3]}, mesh);
}

TP_mesh
double geom_tetra_volume(const key_vector &node_keys, mesh_type &mesh) {
    /// basic shape
    /// top point - 0
    /// bottom tria - 1, 2, 3
    auto &n0 = mesh.get_node(node_keys[0]),
            &n1 = mesh.get_node(node_keys[1]),
            &n2 = mesh.get_node(node_keys[2]),
            &n3 = mesh.get_node(node_keys[3]);
    Vec3D v01 = n1.position - n0.position,
            v02 = n2.position - n0.position,
            v03 = n3.position - n0.position;
    return fabs(v01 * (v02 ^ v03)) / 6.0;
}

TP_mesh
double geom_prism_volume(const key_vector &node_keys, mesh_type &mesh) {
    /// top point - 4
    /// bottom quad - 0, 1, 2, 3
    return geom_tetra_volume<mesh_type>({node_keys[0], node_keys[1], node_keys[3], node_keys[4]}, mesh)
           + geom_tetra_volume<mesh_type>({node_keys[2], node_keys[1], node_keys[3], node_keys[4]}, mesh);
}

TP_mesh
double geom_pyram_volume(const key_vector &node_keys, mesh_type &mesh) {
    /// top - tria 0, 1, 2
    /// bottom - tria 3, 4, 5
    return geom_prism_volume<mesh_type>(
            {node_keys[0], node_keys[2], node_keys[4], node_keys[3], node_keys[2]}, mesh)
           + geom_tetra_volume<mesh_type>({node_keys[2], node_keys[3], node_keys[4], node_keys[5]}, mesh);
}

TP_mesh
double geom_hexah_volume(const key_vector &node_keys, mesh_type &mesh) {
    return geom_pyram_volume<mesh_type>({node_keys[0], node_keys[1], node_keys[3],
                                                   node_keys[4], node_keys[5], node_keys[7]}, mesh)
           + geom_pyram_volume<mesh_type>({node_keys[2], node_keys[1], node_keys[3],
                                                     node_keys[6], node_keys[5], node_keys[7]}, mesh);
}
