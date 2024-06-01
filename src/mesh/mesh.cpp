#include "mesh/mesh.h"


MESO::Mesh::Zone MESO::Mesh::load_gambit(const std::string &file_path) {
    MESO::Mesh::Reader::Gambit reader(file_path);
    MESO::Mesh::Zone mesh = reader.read();
    mesh.build_geom();
    return mesh;
}
