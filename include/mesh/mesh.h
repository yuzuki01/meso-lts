#ifndef MESO_MESH_H
#define MESO_MESH_H

#include "core/core.h"

/// UDF
#include "mesh/object.h"
#include "mesh/geom.h"
#include "mesh/reader.h"

namespace MESO::Mesh {
    MESO::Mesh::Zone load_gambit(const std::string &file_path);
}

#endif //MESO_MESH_H
