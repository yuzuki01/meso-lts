#include "mesh.h"


using namespace DPM;

TP_func bool Particle<int>::is_in_cell(const int &cell_key) {
    return is_in_cell(mesh.get_cell(cell_key));
}

TP_func bool Particle<std::string>::is_in_cell(const std::string &cell_key) {
    return is_in_cell(mesh.get_cell(cell_key));
}
