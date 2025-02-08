#include "mesh/mesh.h"

namespaceMesoMesh


Cell::Cell(const Mesh::GeomMesh &owner,
           const ObjectId &id,
           const ObjectType &geomType,
           const List<ObjectId> &nodes)
        : NodeBasedObject(owner, id, geomType, nodes),
          V_(Geometric::getVolume(mesh_, nodes_, GT_)) {
    // reserve
    neighbors_.reserve(Geometric::faceNum(GT_));
    // resize
    faces_.resize(Geometric::faceNum(GT_), -1);
}

const Scalar &Cell::V() const {
    return V_;
}

void Cell::setFace(const ObjectId &idFace,
                   const ObjectId &idOnOwner) {
    faces_[idOnOwner] = idFace;
}

const List<ObjectId> &Cell::faces() const {
    return faces_;
}

void Cell::setNeighbor(const Label &nei) {
    if (std::find(neighbors_.begin(), neighbors_.end(), nei) == neighbors_.end()) {
        neighbors_.push_back(nei);
    }
}

const List<ObjectId> &Cell::neighbors() const {
    return neighbors_;
}

void Cell::setPart(const Label &rank, const Label &idOnPart) {
    rank_ = rank;
    idOnPartition_ = idOnPart;
}

const Label &Cell::rank() const {
    return rank_;
}

const Label &Cell::idOnPartition() const {
    return idOnPartition_;
}
