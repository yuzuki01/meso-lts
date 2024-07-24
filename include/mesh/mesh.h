#ifndef MESO_MESH_H
#define MESO_MESH_H

#include "core/core.h"

/// METIS
#include "mesh/metis.h"

/// UDF
#include "mesh/object.h"
#include "mesh/geom.h"
#include "mesh/reader.h"

namespace MESO::fvmMesh {
    MESO::fvmMesh::Mesh load_gambit(const std::string &file_path, double mesh_scale=1.0);
}

#endif //MESO_MESH_H
