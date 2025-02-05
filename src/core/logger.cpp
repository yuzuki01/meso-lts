#include "core/core.h"

using namespace MESO;

Logger MESO::logger;


bool Printer::isOutputAvailable() {
    return (root_ptr->level <= output_level) and (MPI::rank == MPI::mainRank);
}

Printer &Printer::operator<<(const Vector &value) {
    if (isOutputAvailable()) {
        std::cout << color_prefix << value.str();
    }
    return *this;
}
