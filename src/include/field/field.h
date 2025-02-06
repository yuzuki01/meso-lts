#ifndef MESO_FIELD_H
#define MESO_FIELD_H

#include "core/core.h"
#include "mesh/mesh.h"


namespaceMesoMesh


#define FieldTemplate template<typename ValueType, Label PatchType>

namespace MESO {

    template<typename ValueType, Label PatchType>
    class BasicField;

    enum {
        VolFlag, SurfFlag
    };

    template<typename ValueType>
    using volField = BasicField<ValueType, VolFlag>;
    using volScalarField = volField<Scalar>;
    using volVectorField = volField<Vector>;

    template<typename ValueType>
    using surfField = BasicField<ValueType, SurfFlag>;
    using surfScalarField = surfField<Scalar>;
    using surfVectorField = surfField<Vector>;

    FieldTemplate
    class BasicField {
    public:
        const Label flag_ = PatchType;
    protected:
        const fvMesh &mesh_;
        const Time &time_;
        List<ObjectId> index_;
        List<ValueType> values_;

    public:
        /// Constructors
        explicit BasicField(const fvMesh &mesh);

        BasicField(const fvMesh &mesh, const List<ValueType> &values);

        ~BasicField() = default;

        /// Interfaces
        List<ValueType> &values();

        [[nodiscard]] const List<ValueType> &values() const;

        [[nodiscard]] const fvMesh &mesh() const;

        [[nodiscard]] const Time &time() const;

        ValueType &operator[](const Label &index);

        const ValueType &operator[](const Label &index) const;

        /// Math

        BasicField &operator=(const BasicField &other);

        BasicField operator+(const BasicField &other);

        BasicField operator-(const BasicField &other);
    };
}

/// Math
template<typename T1, typename T2, Label PatchType>
void isOperationAvailable(const BasicField<T1, PatchType> &_x, const BasicField<T2, PatchType> &_y);

FieldTemplate
BasicField<ValueType, PatchType> operator+(const ValueType &_x, const BasicField<ValueType, PatchType> &_y);

FieldTemplate
BasicField<ValueType, PatchType> operator-(const ValueType &_x, const BasicField<ValueType, PatchType> &_y);

/// Constructors
#include "field/fieldConstructors.h"
/// Interfaces
#include "field/fieldInterfaces.h"
/// Math
#include "field/fieldMath.h"

#endif //MESO_FIELD_H
