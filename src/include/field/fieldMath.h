#ifndef MESO_FIELDMATH_H
#define MESO_FIELDMATH_H

#include "field.h"

template<typename T, Label PatchType>
void isOperationAvailable(const BasicField<T, PatchType> &_x, const BasicField<T, PatchType> &_y) {
    if (_x.mesh().name() == _y.mesh().name()) return;
    logger.error << "BasicField caught unexpected operation" << std::endl;
    FATAL_ERROR_THROW;
}

template<typename T1, typename T2, Label PatchType>
void isOperationAvailable(const BasicField<T1, PatchType> &_x, const BasicField<T2, PatchType> &_y) {
    if (_x.mesh().name() == _y.mesh().name()) return;
    logger.error << "BasicField caught unexpected operation" << std::endl;
    FATAL_ERROR_THROW;
}

FieldTemplate
BasicField<ValueType, PatchType> &
BasicField<ValueType, PatchType>::operator=(const ValueType &x) {
    forAll(values_, i) {
        values_[i] = x;
    }
    return *this;
}

FieldTemplate
BasicField<ValueType, PatchType> &
BasicField<ValueType, PatchType>::operator=(const BasicField<ValueType, PatchType> &other) {
    isOperationAvailable<ValueType, PatchType>(*this, other);
    values_ = other.values_;
    return *this;
}

FieldTemplate
BasicField<ValueType, PatchType> &
BasicField<ValueType, PatchType>::operator+=(const ValueType &x) {
    forAll(values_, i) {
        values_[i] += x;
    }
    return *this;
}

FieldTemplate
BasicField<ValueType, PatchType> &
BasicField<ValueType, PatchType>::operator+=(const BasicField<ValueType, PatchType> &other) {
    isOperationAvailable<ValueType, PatchType>(*this, other);
    forAll(values_, i) {
        values_[i] += other.values()[i];
    }
    return *this;
}

FieldTemplate
BasicField<ValueType, PatchType>
BasicField<ValueType, PatchType>::operator+(const BasicField<ValueType, PatchType> &other) {
    List<ValueType> result(values_.size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = values_[i] + other.values()[i];
    }
    return {mesh_, result};
}


FieldTemplate
BasicField<ValueType, PatchType>
BasicField<ValueType, PatchType>::operator-(const BasicField<ValueType, PatchType> &other) {
    List<ValueType> result(values_.size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = values_[i] - other.values()[i];
    }
    return {mesh_, result};
}

FieldTemplate
BasicField<ValueType, PatchType>
BasicField<ValueType, PatchType>::operator-() {
    List<ValueType> result(values_.size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = -values_[i];
    }
    return {mesh_, result};
}

FieldTemplate
BasicField<ValueType, PatchType> operator+(const ValueType &_x, const BasicField<ValueType, PatchType> &_y) {
    List<ValueType> result(_y.values().size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = _x + _y.values()[i];
    }
    return {_y.mesh(), result};
}


FieldTemplate
BasicField<ValueType, PatchType> operator-(const ValueType &_x, const BasicField<ValueType, PatchType> &_y) {
    List<ValueType> result(_y.values().size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = _x - _y.values()[i];
    }
    return {_y.mesh(), result};
}

PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const Scalar &_x, const BasicField<Scalar, PatchType> &_y) {
    List<Scalar> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Vector, PatchType> operator*(const Vector &_x, const BasicField<Scalar, PatchType> &_y) {
    List<Vector> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const Vector &_x, const BasicField<Vector, PatchType> &_y) {
    List<Scalar> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const BasicField<Scalar, PatchType> &_y, const Scalar &_x) {
    List<Scalar> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Vector, PatchType> operator*(const BasicField<Scalar, PatchType> &_y, const Vector &_x) {
    List<Vector> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const BasicField<Vector, PatchType> &_y, const Vector &_x) {
    List<Scalar> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x * _y.values()[i];
    }
    return {_y.mesh(), result};
}


PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const BasicField<Scalar, PatchType> &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Scalar> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] * _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator*(const BasicField<Scalar, PatchType> &_x,
                                        const BasicField<Vector, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] * _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator*(const BasicField<Vector, PatchType> &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] * _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Scalar, PatchType> operator*(const BasicField<Vector, PatchType> &_x,
                                        const BasicField<Vector, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Scalar> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] * _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator^(const BasicField<Vector, PatchType> &_x,
                                        const BasicField<Vector, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] ^ _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Scalar, PatchType> operator/(const BasicField<Scalar, PatchType> &_x,
                                        const Scalar &_y) {
    List<Scalar> result(_x.values().size());
    Scalar _s = 1.0 / _y;
    forAll(result, i) {
        result[i] = _x.values()[i] * _s;
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator/(const BasicField<Vector, PatchType> &_x,
                                        const Scalar &_y) {
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] / _y;
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Scalar, PatchType> operator/(const Scalar &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    List<Scalar> result(_y.values().size());
    Scalar _s = 1.0 / _x;
    forAll(result, i) {
        result[i] = _y.values()[i] * _s;
    }
    return {_y.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator/(const Vector &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    List<Vector> result(_y.values().size());
    forAll(result, i) {
        result[i] = _x / _y.values()[i];
    }
    return {_y.mesh(), result};
}

PatchTypeTemplate
BasicField<Scalar, PatchType> operator/(const BasicField<Scalar, PatchType> &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] / _y.values()[i];
    }
    return {_x.mesh(), result};
}

PatchTypeTemplate
BasicField<Vector, PatchType> operator/(const BasicField<Vector, PatchType> &_x,
                                        const BasicField<Scalar, PatchType> &_y) {
    isOperationAvailable(_x, _y);
    List<Vector> result(_x.values().size());
    forAll(result, i) {
        result[i] = _x.values()[i] / _y.values()[i];
    }
    return {_x.mesh(), result};
}

#endif //MESO_FIELDMATH_H
