#include <utility>

#include "mesh/mesh.h"


namespaceMesoMesh

Patch::Patch(const Mesh::GeomMesh &owner,
             String name)
           : mesh_(owner), name_(std::move(name)) {

}

void Patch::append(const ObjectId &objectId) {
    group_.push_back(objectId);
}

void Patch::fit() {
    group_.shrink_to_fit();
}

Label Patch::size() const {
    return static_cast<Label>(group_.size());
}

const String &Patch::name() const {
    return name_;
}

const GeomMesh &Patch::mesh() const {
    return mesh_;
}

ObjectId &Patch::operator[](const Label &index) {
    return group_[index];
}

const ObjectId &Patch::operator[](const Label &index) const {
    return group_[index];
}

List<ObjectId> &Patch::group() {
    return group_;
}

const List<ObjectId> &Patch::group() const {
    return group_;
}
