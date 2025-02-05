#include <utility>

#include "mesh/mesh.h"


namespaceMesoMesh

Patch::Patch(const MESO::Mesh::GeomMesh &owner,
             MESO::String name)
           : mesh_(owner), name(std::move(name)) {

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
