#include <utility>

#include "mesh/mesh.h"


namespaceMesoMesh

Patch::Patch(const Mesh::GeomMesh &owner,
             String name)
           : mesh_(owner), name_(std::move(name)) {

}

const String &Patch::name() const {
    return name_;
}

const GeomMesh &Patch::mesh() const {
    return mesh_;
}

List<ObjectId> &Patch::group() {
    return group_;
}

const List<ObjectId> &Patch::group() const {
    return group_;
}
