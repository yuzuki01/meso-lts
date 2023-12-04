/**
 * included by mesh.h
 */

#define MESH_GAUSS_HERMIT "GaussHermit"
#define MESH_NEWTON_COTES "NewtonCotes"

namespace GENERATOR {
    MESH::StaticMesh gauss_hermit(int gauss_point, int dimension, double RT);
    MESH::StaticMesh newton_cotes(int n, int mount, int dimension, double scale, double RT);
}
