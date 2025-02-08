#ifndef MESO_FIELDMATHFUNCTIONS_H
#define MESO_FIELDMATHFUNCTIONS_H

namespace MESO {

    PatchTypeTemplate
    BasicField<Scalar, PatchType> sqrt(const BasicField<Scalar, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = std::sqrt(_x[i]);
        }
        return BasicField(_x.mesh(), result);
    }

    PatchTypeTemplate
    BasicField<Scalar, PatchType> exp(const BasicField<Scalar, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = std::exp(_x[i]);
        }
        return BasicField(_x.mesh(), result);
    }

    PatchTypeTemplate
    BasicField<Scalar, PatchType> mag(const BasicField<Vector, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = _x[i].magnitude();
        }
        return BasicField(_x.mesh(), result);
    }
}

#endif //MESO_FIELDMATHFUNCTIONS_H
