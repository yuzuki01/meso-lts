#include "mesh/mesh.h"

namespaceMesoMesh

Node::Node(
        const Mesh::GeomMesh &owner,
        const ObjectId &id,
        const Vector &coordinate)
        : BasicMeshObject(owner, id, Geometric::Node),
        C_(coordinate) {

}

const Coordinate &Node::C() const {
    return C_;
}
