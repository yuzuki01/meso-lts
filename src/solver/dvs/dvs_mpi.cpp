#include "solver/solver.h"

using namespace MESO;

MPI::MPI_TaskObject MPI::DVS_partition(Mesh::Mesh &mesh) {
    auto task_list = MPI::get_task_distribution(mesh.NCELL);
    mesh.cell_names.resize(MPI::node_num);
    mesh.cell_groups.resize(MPI::node_num);
    for (int i = 0; i < MPI::node_num; ++i) {
        auto &task = task_list[i];
        auto &group = mesh.cell_groups[i];
        group.resize(task.size);
        mesh.cell_names[i] = "node-" + std::to_string(MPI::rank);
        for (int j = 0; j < task.size; ++j) {
            group[j] = task.start + j;
        }
    }
    return task_list[MPI::rank];
}
