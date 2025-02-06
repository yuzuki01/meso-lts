#include "mesh/mesh.h"

namespaceMesoMesh


Cell::Cell(const Mesh::GeomMesh &owner,
           const ObjectId &id,
           const ObjectType &geomType,
           const List<ObjectId> &nodes)
        : NodeBasedObject(owner, id, geomType, nodes),
          V_(Geometric::getVolume(mesh_, nodes_, GT_)) {
    // resize
    faces_.resize(Geometric::faceNum(GT_), -1);
}

void Cell::setFace(const ObjectId &idFace,
                   const ObjectId &idOnOwner) {
    faces_[idOnOwner] = idFace;
}

const Scalar &Cell::V() const {
    return V_;
}

const List<ObjectId> &Cell::faces() const {
    return faces_;
}
