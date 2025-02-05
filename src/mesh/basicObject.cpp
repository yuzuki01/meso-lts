#include "mesh/mesh.h"

namespaceMesoMesh

BasicMeshObject::BasicMeshObject(
        const GeomMesh &owner,
        const ObjectId &id,
        const ObjectType& GT)
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



