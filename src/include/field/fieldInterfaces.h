#ifndef MESO_FIELDINTERFACES_H
#define MESO_FIELDINTERFACES_H

FieldTemplate
List<ValueType> &BasicField<ValueType, PatchType>::values() {
    return values_;
}

FieldTemplate
[[nodiscard]] const List<ObjectId> &BasicField<ValueType, PatchType>::index() const {
    return index_;
}

FieldTemplate
[[nodiscard]] const List<ValueType> &BasicField<ValueType, PatchType>::values() const {
    return values_;
}

FieldTemplate
[[nodiscard]] const Label BasicField<ValueType, PatchType>::size() const {
    return values_.size();
}

FieldTemplate
[[nodiscard]] const fvMesh &BasicField<ValueType, PatchType>::mesh() const {
    return mesh_;
}

FieldTemplate
[[nodiscard]] const Time &BasicField<ValueType, PatchType>::time() const {
    return time_;
}

FieldTemplate
ValueType &BasicField<ValueType, PatchType>::operator[](const Label &index) {
    return values_[index];
}

FieldTemplate
const ValueType &BasicField<ValueType, PatchType>::operator[](const Label &index) const {
    return values_[index];
}

#endif //MESO_FIELDINTERFACES_H
