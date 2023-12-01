#include "mesh.h"


using namespace DPM;

TP_func bool Particle<MESH::StaticMesh>::is_in_cell(int cell_key) {
    return is_in_cell(mesh.get_cell(cell_key));
}

TP_func bool Particle<MESH::MapMesh>::is_in_cell(int cell_key) {
    return is_in_cell(mesh.get_cell(cell_key));
}
