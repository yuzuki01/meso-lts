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

#endif //MESO_FIELDINTERFACES_H
