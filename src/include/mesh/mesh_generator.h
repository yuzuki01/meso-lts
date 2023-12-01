/**
 * included by mesh.h
 */

#define MESH_GAUSS_HERMIT "GaussHermit"

namespace GENERATOR {
    MESH::StaticMesh gauss_hermit(int gauss_point, int dimension, double RT);
}
