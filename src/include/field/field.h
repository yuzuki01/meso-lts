#ifndef MESO_FIELD_H
#define MESO_FIELD_H

#include "core/core.h"
#include "mesh/mesh.h"


namespaceMesoMesh

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

    template<typename ValueType, Label PatchType>
    class BasicField {
    protected:
        const fvMesh &mesh_;
        const Time &time_;
        List<ValueType> data_;

    public:
        explicit BasicField(const fvMesh &mesh);
    };
}

#endif //MESO_FIELD_H
