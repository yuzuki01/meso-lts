#ifndef MESO_FIELD_H
#define MESO_FIELD_H

#include "core/core.h"
#include "mesh/mesh.h"


namespaceMesoMesh


#define FieldTemplate template<typename ValueType, Label PatchType>
#define PatchTypeTemplate template<Label PatchType>


namespace MESO {

    FieldTemplate
    class BasicField {
    public:
        const Label flag_ = PatchType;
    protected:
        const fvMesh &mesh_;
        const Time &time_;
        List<ObjectId> index_;
        List<ValueType> values_;

        /// Constructors
    public:
        explicit BasicField(const fvMesh &mesh);

        BasicField(const fvMesh &mesh,
                   const List<ValueType> &values);

        ~BasicField() = default;

        /// Interfaces
    public:
        List<ValueType> &values();

        [[nodiscard]] const List<ObjectId> &index() const;

        [[nodiscard]] const List<ValueType> &values() const;

        [[nodiscard]] const Label size() const;

        [[nodiscard]] const fvMesh &mesh() const;

        [[nodiscard]] const Time &time() const;

        BasicField<Scalar, PatchType> component(const Label &index) const;

        /// Math
    public:
        BasicField &operator=(const BasicField &other);

        BasicField operator+(const BasicField &other);

        BasicField operator-(const BasicField &other);

        /// File IO
        void output(const String &fieldName);
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
#include "field/fieldMathFunctions.h"
/// fvMesh
#include "field/fvMeshField.h"

#endif //MESO_FIELD_H
