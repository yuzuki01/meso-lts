#include "field/field.h"


namespaceMESO


template<>
BasicField<Scalar, VolFlag> BasicField<Vector, VolFlag>::component(const MESO::Label &index) const {
    List<Scalar> result(values_.size());
    switch (index) {
        case X_COMPONENT:
        case Y_COMPONENT:
        case Z_COMPONENT:
            forAll(values_, i)
            {
                result[i] = values_[i][index];
            }
            return {mesh_, result};
        default:
            logger.error << "VectorField::component() caught index out of range" << std::endl;
            FATAL_ERROR_THROW;
    }
}


template<>
BasicField<Scalar, SurfFlag> BasicField<Vector, SurfFlag>::component(const MESO::Label &index) const {
    List<Scalar> result(values_.size());
    switch (index) {
        case X_COMPONENT:
        case Y_COMPONENT:
        case Z_COMPONENT:
            forAll(values_, i)
            {
                result[i] = values_[i][index];
            }
            return {mesh_, result};
        default:
            logger.error << "VectorField::component() caught index out of range" << std::endl;
            FATAL_ERROR_THROW;
    }
}
