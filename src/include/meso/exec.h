#ifndef MESO_MPI_EXEC_H
#define MESO_MPI_EXEC_H

namespace MESO {
    void help();

    int handle_mesh(const std::string &mesh_file, double mesh_scale=1.0);
}

#endif //MESO_MPI_EXEC_H
