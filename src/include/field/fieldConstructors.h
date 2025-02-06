#ifndef MESO_FIELDCONSTRUCTOR_H
#define MESO_FIELDCONSTRUCTOR_H


template<>
// volField
explicit BasicField<ValueType, VolFlag>(const fvMesh &mesh)
        : mesh_(mesh), time_(mesh.time()) {
    const auto &partitionPatch = mesh_.partition()[MPI::rank];
    const auto size = static_cast<Label>(partitionPatch.size());
    values_.resize(size, ValueType(0));
    for (int i = 0; i < size; ++i) {
        index_.push_back(partitionPatch[i]);
    }
}

template<>
// surfField
explicit BasicField<ValueType, SurfFlag>(const fvMesh &mesh)
        : mesh_(mesh), time_(mesh_.time()) {
    const auto &partitionPatch = mesh_.partition()[MPI::rank];
    const auto size = static_cast<Label>(partitionPatch.size());
    for (int ci = 0; ci < size; ++ci) {
        const auto &cell = mesh_.cell(ci);
        for (const auto &fi: cell.faces()) {
            const auto &face = mesh_.face(fi);
            if (
                // face is owned by cell
                    (face.owner() == ci)
                    and
                    // faceId does not exist
                    (std::find(index_.begin(), index_.end(), fi) == index_.end())
                    ) {
                index_.push_back(fi);
            }
        }
    }
    index_.shrink_to_fit();
    values_.resize(index_.size(), ValueType(0));
}

#endif //MESO_FIELDCONSTRUCTOR_H
