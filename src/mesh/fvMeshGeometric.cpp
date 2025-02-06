#include "mesh/mesh.h"

namespaceMesoMesh

namespace Area {
    Scalar edge(const GeomMesh &mesh, const List<ObjectId> &nodes);

    Scalar tria(const GeomMesh &mesh, const List<ObjectId> &nodes);

    Scalar quad(const GeomMesh &mesh, const List<ObjectId> &nodes);
}

namespace Volume {
    using Area::tria;
    using Area::quad;

    Scalar brick(const GeomMesh &mesh, const List<ObjectId> &nodes);

    Scalar wedge(const GeomMesh &mesh, const List<ObjectId> &nodes);

    Scalar tetra(const GeomMesh &mesh, const List<ObjectId> &nodes);

    Scalar pyram(const GeomMesh &mesh, const List<ObjectId> &nodes);
}

namespace Spatial {
    inline Scalar vectorAngleXOY(const Vector &_vec1, const Vector &_vec2);
}

/**
 * =================================================
 * -------------------- Global ---------------------
 * =================================================
 **/

Label Geometric::nodeNum(const ObjectType &geomType) {
    switch (geomType) {
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
            logger.error << "Geometric::nodeNum() caught unexpected value" << std::endl;
            FATAL_ERROR_THROW;
    }
}

Label Geometric::faceNum(const ObjectType &geomType) {
    switch (geomType) {
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
            logger.error << "Geometric::faceNum() caught unexpected value" << std::endl;
            FATAL_ERROR_THROW;
    }
}

Vector Geometric::getCoordinate(const GeomMesh &mesh,
                                const List<ObjectId> &nodes) {
    Vector C(0, 0, 0);
    for (const auto &ni: nodes) {
        auto &node = mesh.nodes()[ni];
        C += node.C();
    }
    return C / static_cast<Scalar>(nodes.size());
}

Scalar Geometric::getArea(const GeomMesh &mesh,
                          const List<ObjectId> &nodes,
                          const ObjectType &GT) {
    switch (GT) {
        case Edge:
            return Area::edge(mesh, nodes);
        case Tria:
            return Area::tria(mesh, nodes);
        case Quad:
            return Area::quad(mesh, nodes);
        default:
            logger.error << "Geometric::getArea() caught unexpected GeomType = "
                         << GT << std::endl;
            FATAL_ERROR_THROW;
    }
}

Vector Geometric::getSurfaceVector(const GeomMesh &mesh,
                                   const List<ObjectId> &nodes,
                                   const ObjectType &GT,
                                   const ObjectId &owner,
                                   const Vector &faceC,
                                   const Scalar &Sf) {
    auto &own = mesh.cells()[owner];
    List<const Mesh::Node *> pNodes;
    for (const auto &ni: nodes) {
        auto &node = mesh.node(ni);
        pNodes.push_back(&node);
    }
    Vector vcf = faceC - mesh.cell(owner).C();
    Vector vfn1 = (pNodes[0]->C() - pNodes[1]->C()).normalize();
    Vector Snv;
    if (GT == Edge) {
        // 2D
        auto aglX = Spatial::vectorAngleXOY(vfn1, {1.0, 0.0, 0.0});
        auto aglY = Spatial::vectorAngleXOY(vfn1, {0.0, 1.0, 0.0});
        if (aglY <= M_PI_2) {
            aglX = M_PI_2 + aglX;
        } else {
            aglX = M_PI_2 - aglX;
        }
        Snv = Vector(cos(aglX), sin(aglX), 0.0).normalize() * Sf;
    } else {
        // 3D
        Vector vfn2 = (pNodes[1]->C() - pNodes[2]->C()).normalize();
        Snv = (vfn1 ^ vfn2).normalize() * Sf;
    }
    return ((Snv * vcf) >= 0) ? Snv : -Snv;
}

Scalar Geometric::getVolume(const GeomMesh &mesh,
                            const List<ObjectId> &nodes,
                            const ObjectType &GT) {
    switch (GT) {
        case Tria:
            return Volume::tria(mesh, nodes);
        case Quad:
            return Volume::quad(mesh, nodes);
        case Brick:
            return Volume::brick(mesh, nodes);
        case Wedge:
            return Volume::wedge(mesh, nodes);
        case Tetra:
            return Volume::tetra(mesh, nodes);
        case Pyram:
            return Volume::pyram(mesh, nodes);
        default:
            logger.error << "Geometric::getVolume() caught unexpected GeomType = "
                         << GT << std::endl;
            FATAL_ERROR_THROW;
    }
}

/**
 * =================================================
 * --------------------- Local ---------------------
 * =================================================
 **/

Scalar Spatial::vectorAngleXOY(const Vector &_vec1, const Vector &_vec2) {
    Scalar m1 = _vec1.magnitude(), m2 = _vec2.magnitude();
    return acos((_vec1 * _vec2) / (m1 * m2));
}

inline Vector constructVector(const Node &n1, const Node &n2) {
    return n2.C() - n1.C();
}

// Area

Scalar Area::edge(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    return constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[1])).magnitude();
}

Scalar Area::tria(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    Vector p1 = constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[1]));
    Vector p2 = constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[2]));
    return 0.5 * (p1 ^ p2).magnitude();
}

Scalar Area::quad(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    return tria(mesh, {nodes[0], nodes[1], nodes[2]})
           + tria(mesh, {nodes[0], nodes[3], nodes[2]});
}

Scalar Volume::brick(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    return wedge(mesh, {nodes[0], nodes[1], nodes[5], nodes[2], nodes[3], nodes[7]})
           + wedge(mesh, {nodes[0], nodes[4], nodes[5], nodes[2], nodes[6], nodes[7]});
}

Scalar Volume::wedge(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    return pyram(mesh, {nodes[0], nodes[2], nodes[5], nodes[3], nodes[1]})
           + tetra(mesh, {nodes[1], nodes[3], nodes[4], nodes[5]});
}

Scalar Volume::pyram(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    return tetra(mesh, {nodes[0], nodes[1], nodes[2], nodes[4]})
           + tetra(mesh, {nodes[1], nodes[2], nodes[3], nodes[4]});
}

Scalar Volume::tetra(const GeomMesh &mesh, const List<ObjectId> &nodes) {
    Vector p1 = constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[1]));
    Vector p2 = constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[2]));
    Vector p3 = constructVector(
            mesh.node(nodes[0]),
            mesh.node(nodes[3]));
    return std::abs(p1 * (p2 ^ p3)) / 6.0;
}
