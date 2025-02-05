#include "mesh/mesh.h"

namespaceMesoMesh

BasicMeshObject::BasicMeshObject(
        const GeomMesh &owner,
        const ObjectId &id,
        const ObjectType &GT)
        : mesh_(owner), id_(id), GT_(GT) {}

const GeomMesh &BasicMeshObject::mesh() const {
    return mesh_;
}

const ObjectId &BasicMeshObject::id() const {
    return id_;
}

const ObjectType &BasicMeshObject::geomType() const {
    return GT_;
}


NodeBasedObject::NodeBasedObject(
        const MESO::Mesh::GeomMesh &owner,
        const MESO::ObjectId &id,
        const MESO::ObjectType &GT,
        const List<MESO::ObjectId> &nodes)
        : BasicMeshObject(owner, id, GT),
          nodes_(nodes),
          C_(Geometric::getCoordinate(mesh_, nodes)) {

}

const List<ObjectId> &NodeBasedObject::nodes() const {
    return nodes_;
}

const Vector &NodeBasedObject::C() const {
    return C_;
}
