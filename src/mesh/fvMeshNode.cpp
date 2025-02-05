#include "mesh/mesh.h"

namespaceMesoMesh

Node::Node(
        const MESO::Mesh::GeomMesh &owner,
        const MESO::ObjectId &id,
        const MESO::Vector &coordinate)
        : BasicMeshObject(owner, id, Geometric::Node),
        C_(coordinate) {

}

const Vector &Node::C() const {
    return C_;
}
