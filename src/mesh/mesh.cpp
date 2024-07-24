#include "mesh/mesh.h"


MESO::fvmMesh::Mesh MESO::fvmMesh::load_gambit(const std::string &file_path, double mesh_scale) {
    logger.info << "Load gambit mesh: " << file_path << std::endl;
    MESO::fvmMesh::Reader::Gambit reader(file_path);
    MESO::fvmMesh::Mesh mesh = reader.parse_mesh_strings();
    mesh.build_geom(mesh_scale);
    logger.info << "Built mesh(" << file_path << ") topology." << std::endl;
    return mesh;
}
