#include "mesh/mesh.h"


MESO::Mesh::Mesh MESO::Mesh::load_gambit(const std::string &file_path) {
    logger.info << "Load gambit mesh: " << file_path << std::endl;
    MESO::Mesh::Reader::Gambit reader(file_path);
    MESO::Mesh::Mesh mesh = reader.read();
    mesh.build_geom();
    logger.info << "Built mesh(" << file_path << ") topology." << std::endl;
    return mesh;
}
