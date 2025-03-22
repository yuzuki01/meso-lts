#include "mesh/mesh.h"


using namespace MESO;

Vector edge_vector(int start, int end, const List<fvmMesh::Node> &node_list);

double vector_angle_2d(const Vector &_vec1, const Vector &_vec2);

namespace Area {
    double edge(const List<fvmMesh::Node> &node_list);

    double tria(const List<fvmMesh::Node> &node_list);

    double quad(const List<fvmMesh::Node> &node_list);
}

namespace Volume {
    using Area::tria;
    using Area::quad;

    double brick(const List<fvmMesh::Node> &node_list);

    double wedge(const List<fvmMesh::Node> &node_list);

    double tetra(const List<fvmMesh::Node> &node_list);

    double pyram(const List<fvmMesh::Node> &node_list);
}

int Geom::face_num(ObjectType geom_type) {
    switch (geom_type) {
        case Tria:
            return 3;
        case Quad:
        case Tetra:
            return 4;
        case Wedge:
        case Pyram:
            return 5;
        case Brick:
            return 6;
        default:
            throw std::invalid_argument("<MESO::Geom::face_num> unsupported geom_type");
    }
}

int Geom::node_num(MESO::ObjectType geom_type) {
    switch (geom_type) {
        case Edge:
            return 2;
        case Tria:
            return 3;
        case Quad:
        case Tetra:
            return 4;
        case Pyram:
            return 5;
        case Wedge:
            return 6;
        case Brick:
            return 8;
        default:
            throw std::invalid_argument("<MESO::Geom::node_num> unsupported geom_type");
    }
}

Vector Geom::calculate_position(const MESO::List<fvmMesh::Node> &node_list) {
    Vector pos(0.0, 0.0, 0.0);
    const int count = int(node_list.size());
    for (auto &it: node_list) {
        pos += it.position;
    }
    return pos * (1.0 / count);
}

Scalar Geom::calculate_area(ObjectType geom_type, const List<fvmMesh::Node> &node_list) {
    switch (geom_type) {
        case Edge:
            return Area::edge(node_list);
        case Tria:
            return Area::tria(node_list);
        case Quad:
            return Area::quad(node_list);
        default:
            throw std::invalid_argument("<MESO::Geom::calculate_area> unsupported geom_type");
    }
}

Scalar Geom::calculate_volume(MESO::ObjectType geom_type, const MESO::List<fvmMesh::Node> &node_list) {
    switch (geom_type) {
        case Tria:
            return Volume::tria(node_list);
        case Quad:
            return Volume::quad(node_list);
        case Brick:
            return Volume::brick(node_list);
        case Wedge:
            return Volume::wedge(node_list);
        case Tetra:
            return Volume::tetra(node_list);
        case Pyram:
            return Volume::pyram(node_list);
        default:
            throw std::invalid_argument("<Geom::calculate_volume> unsupported geom_type");
    }
}

Vector Geom::calculate_normal_vector(ObjectType face_geom_type, MESO::Vector &cell_position,
                                     MESO::Vector &face_position, const MESO::List<fvmMesh::Node> &node_list) {
    Vector vcf = face_position - cell_position,
            vfn = edge_vector(0, 1, node_list).normalize();
    if (face_geom_type == Edge) {
        /// 2D
        /// get angle
        double angle_x = vector_angle_2d(vfn, {1.0, 0.0, 0.0});
        double angle_y = vector_angle_2d(vfn, {0.0, 1.0, 0.0});
        if (angle_y <= M_PI_2) {
            angle_x = M_PI_2 + angle_x;
        } else {
            angle_x = M_PI_2 - angle_x;
        }
        /// get normal vector
        Vector nv = Vector(cos(angle_x), sin(angle_x), 0.0).normalize();
        /// direction
        return (nv * vcf >= 0.0) ? nv : -nv;
    }
    /// 3D
    Vector nv = (edge_vector(1, 0, node_list) ^ edge_vector(1, 2, node_list)).normalize();
    return (nv * vcf >= 0.0) ? nv : -nv;
}

KeyString Geom::generate_key(const MESO::List<ObjectId> &node_list) {
    const int len = int(node_list.size());
    /// find where is min
    int pos = 0;
    int min = node_list[pos];
    for (int i = 0; i < len; ++i) {
        if (node_list[i] < min) {
            pos = i;
            min = node_list[i];
        }
    }
    /// direction
    bool direction;
    int left_pos, right_pos;
    left_pos = (pos == 0) ? (len - 1) : (pos - 1);
    right_pos = (pos == len - 1) ? 0 : (pos + 1);
    direction = node_list[right_pos] < node_list[left_pos];
    /// append key
    int count = 0;
    std::stringstream ss;
    while (count < len) {
        if (count == 0) ss << node_list[pos]; else ss << ">" << node_list[pos];
        pos = direction ? (pos + 1) : (pos - 1);
        /// make valid range
        if (pos < 0) pos = len - 1;
        if (pos > len - 1) pos = 0;
        count++;
    }
    return ss.str();
}


void fvmMesh::Cell::compute_least_square(Mesh &mesh, int dimension) {
    least_square.neighbor_num = int(neighbors.size());
    least_square.dr.resize(least_square.neighbor_num);
    least_square.weight.resize(least_square.neighbor_num);
    double Sxx, Sxy, Sxz, Syy, Syz, Szz;
    Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0;
    for (int i = 0; i < least_square.neighbor_num; i++) {
        auto &neighbor_cell = mesh.cells[neighbors[i]];
        Vector dr = neighbor_cell.position - position;
        double wi = 1.0 / (dr * dr);
        least_square.dr[i] = dr;
        least_square.weight[i] = wi;
        // sum
        Sxx += wi * (dr.x * dr.x);
        Sxy += wi * (dr.x * dr.y);
        Sxz += wi * (dr.x * dr.z);
        Syy += wi * (dr.y * dr.y);
        Syz += wi * (dr.y * dr.z);
        Szz += wi * (dr.z * dr.z);
    }
    least_square.dr.shrink_to_fit();
    least_square.weight.shrink_to_fit();
    if (dimension == 2) {
        double FM = Sxx * Syy - Sxy * Sxy;
        least_square.Cx = {Syy / FM, -Sxy / FM, 0.0};
        least_square.Cy = {-Sxy / FM, Sxx / FM, 0.0};
        least_square.Cz = {0.0, 0.0, 0.0};
    } else {
        double FM = -Sxz * Sxz * Syy + 2.0 * Sxy * Sxz * Syz - Sxx * Syz * Syz - Sxy * Sxy * Szz +
                    Sxx * Syy * Szz;
        double SxyyzNxzyy = Sxy * Syz - Sxz * Syy;
        double SxzyzNxyzz = Sxz * Syz - Sxy * Szz;
        double SyyzzNyzyz = Syy * Szz - Syz * Syz;
        double SxxzzNxzxz = Sxx * Szz - Sxz * Sxz;
        double SxyxzNxxyz = Sxy * Sxz - Sxx * Syz;
        double SxxyyNxyxy = Sxx * Syy - Sxy * Sxy;
        least_square.Cx = {
                SyyzzNyzyz / FM, SxzyzNxyzz / FM, SxyyzNxzyy / FM
        };
        least_square.Cy = {
                SxzyzNxyzz / FM, SxxzzNxzxz / FM, SxyxzNxxyz / FM
        };
        least_square.Cz = {
                SxyyzNxzyy / FM, SxyxzNxxyz / FM, SxxyyNxyxy / FM
        };
    }
}

/// Local
Vector edge_vector(int start, int end, const List<fvmMesh::Node> &node_list) {
    return node_list[end].position - node_list[start].position;
}

double vector_angle_2d(const Vector &_vec1, const Vector &_vec2) {
    double m1 = _vec1.magnitude(), m2 = _vec2.magnitude();
    return acos((_vec1 * _vec2) / (m1 * m2));
}

double Area::edge(const List<fvmMesh::Node> &node_list) {
    return edge_vector(0, 1, node_list).magnitude();
}

double Area::tria(const List<fvmMesh::Node> &node_list) {
    Vector p1 = edge_vector(0, 1, node_list),
            p2 = edge_vector(0, 2, node_list);
    return 0.5 * (p1 ^ p2).magnitude();
}

double Area::quad(const List<fvmMesh::Node> &node_list) {
    return tria({node_list[0], node_list[1], node_list[2]})
           + tria({node_list[0], node_list[3], node_list[2]});
}

double Volume::brick(const MESO::List<fvmMesh::Node> &node_list) {
    return wedge({node_list[0], node_list[1], node_list[5], node_list[2], node_list[3], node_list[7]})
           + wedge({node_list[0], node_list[4], node_list[5], node_list[2], node_list[6], node_list[7]});
}

double Volume::wedge(const MESO::List<fvmMesh::Node> &node_list) {
    return pyram({node_list[0], node_list[2], node_list[5], node_list[3], node_list[1]})
           + tetra({node_list[1], node_list[3], node_list[4], node_list[5]});
}

double Volume::pyram(const MESO::List<fvmMesh::Node> &node_list) {
    return tetra({node_list[0], node_list[1], node_list[2], node_list[4]})
           + tetra({node_list[1], node_list[2], node_list[3], node_list[4]});
}

double Volume::tetra(const MESO::List<fvmMesh::Node> &node_list) {
    Vector p1 = edge_vector(0, 1, node_list),
            p2 = edge_vector(0, 2, node_list),
            p3 = edge_vector(0, 3, node_list);
    return std::abs(p1 * (p2 ^ p3)) / 6.0;
}
