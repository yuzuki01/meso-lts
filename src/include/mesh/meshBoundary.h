#ifndef MESO_MESHBOUNDARY_H
#define MESO_MESHBOUNDARY_H

namespace Boundary {
    enum {
        FluidInterior,
        Processor,
        Wall,
        FarField,
    };

    const Dict<ObjectType> TypeName = {
            {"fluidInterior", FluidInterior},
            {"processor", Processor},
            {"wall", Wall},
            {"farField", FarField},
    };
}

#endif //MESO_MESHBOUNDARY_H
