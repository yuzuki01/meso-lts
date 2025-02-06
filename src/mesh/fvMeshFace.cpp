#include "mesh/mesh.h"

namespaceMesoMesh

Face::Face(const Mesh::GeomMesh &owner,
           const ObjectId &id,
           const ObjectType &geomType,
           const List<ObjectId> &nodes,
           const ObjectId &ci)
        : NodeBasedObject(owner, id, geomType, nodes),
          S_(Geometric::getArea(mesh_, nodes_, GT_)),
          owner_(ci),
          neighbor_(ci),
          Snv_(Geometric::getSurfaceVector(mesh_, nodes_,
                                           GT_, owner_,
                                           C_, S_)) {

}

void Face::setPatch(const ObjectId &patchId) {
    patch_ = patchId;
}

void Face::setNeighbor(const Label &cellId) {
    neighbor_ = cellId;
}

const Scalar &Face::S() const {
    return S_;
}

const ObjectId &Face::owner() const {
    return owner_;
}

const ObjectId &Face::neighbor() const {
    return neighbor_;
}

const ObjectId &Face::patch() const {
    return patch_;
}

const Vector &Face::Snv() const {
    return Snv_;
}
