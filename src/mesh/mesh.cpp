#include "mesh/mesh.h"


MESO::Mesh::Mesh MESO::Mesh::load_gambit(const std::string &file_path) {
    MESO::Mesh::Reader::Gambit reader(file_path);
    MESO::Mesh::Mesh mesh = reader.read();
    mesh.build_geom();
    return mesh;
}
