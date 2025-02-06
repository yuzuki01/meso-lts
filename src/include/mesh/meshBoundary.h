#ifndef MESO_MESHBOUNDARY_H
#define MESO_MESHBOUNDARY_H

namespace Boundary {
    enum {
        FluidInterior,
        Processor,
        Wall,
        FarField,
    };

    const Map<String> TypeName = {
            {FluidInterior, "fluidInterior"},
            {Processor, "processor"},
            {Wall, "wall"},
            {FarField, "farField"},
    };
}

#endif //MESO_MESHBOUNDARY_H
