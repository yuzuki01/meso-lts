#ifndef MESO_FIELDCONSTRUCTOR_H
#define MESO_FIELDCONSTRUCTOR_H

FieldTemplate
BasicField<ValueType, PatchType>::BasicField(const MESO::Mesh::fvMesh &mesh)
        : mesh_(mesh), time_(mesh.time()) {
    switch (PatchType) {
        case VolFlag:   // vol
        {
            const auto &partitionPatch = mesh_.partitionCellPatch()[MPI::rank];
            const auto size = static_cast<Label>(partitionPatch.size());
            for (int i = 0; i < size; ++i) {
                index_.push_back(partitionPatch[i]);
            }
        }
            break;
        case SurfFlag:  // surf
        {
            const auto &partitionPatch = mesh_.partitionFacePatch()[MPI::rank];
            const auto size = static_cast<Label>(partitionPatch.size());
            for (int i = 0; i < size; ++i) {
                index_.push_back(partitionPatch[i]);
            }
        }
            break;
        default:
            logger.error << "BasicField caught unexpected PatchType" << std::endl;
            FATAL_ERROR_THROW;
    }
    values_.resize(index_.size(), ValueType(0));
    index_.shrink_to_fit();
    values_.shrink_to_fit();
}

FieldTemplate
BasicField<ValueType, PatchType>::BasicField(const fvMesh &mesh,
                                             const ValueType &value)
        : BasicField(mesh) {
    forAll(values_, i) {
        values_[i] = value;
    }
}

FieldTemplate
BasicField<ValueType, PatchType>::BasicField(const fvMesh &mesh,
                                             const List<ValueType> &values)
        : BasicField(mesh) {
    values_ = values;
}

#endif //MESO_FIELDCONSTRUCTOR_H
