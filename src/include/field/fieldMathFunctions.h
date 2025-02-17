#ifndef MESO_FIELDMATHFUNCTIONS_H
#define MESO_FIELDMATHFUNCTIONS_H

namespace MESO {

    PatchTypeTemplate
    BasicField<Scalar, PatchType> sqrt(const BasicField<Scalar, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = std::sqrt(_x[i]);
        }
        return BasicField<Scalar, PatchType>(_x.mesh(), _x.index(), result);
    }

    PatchTypeTemplate
    BasicField<Scalar, PatchType> exp(const BasicField<Scalar, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = std::exp(_x.values()[i]);
        }
        return BasicField<Scalar, PatchType>(_x.mesh(), _x.index(), result);
    }

    PatchTypeTemplate
    BasicField<Scalar, PatchType> mag(const BasicField<Vector, PatchType> &_x) {
        List<Scalar> result(_x.values().size());
        forAll(result, i) {
            result[i] = _x[i].magnitude();
        }
        return BasicField<Scalar, PatchType>(_x.mesh(), _x.index(), result);
    }
}

#endif //MESO_FIELDMATHFUNCTIONS_H
