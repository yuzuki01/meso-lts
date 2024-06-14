#include "solver/solver.h"

using namespace MESO;

MESO::MPI::MPI_TaskObject MESO::MPI::DVS_partition(Mesh::Mesh &mesh) {
    auto task_list = MPI::get_task_distribution(mesh.NCELL);
    mesh.cell_names.resize(MPI::process_num);
    mesh.cell_groups.resize(MPI::process_num);
    for (int i = 0; i < MPI::process_num; ++i) {
        auto &task = task_list[i];
        auto &group = mesh.cell_groups[i];
        group.resize(task.size);
        mesh.cell_names[i] = "node-" + std::to_string(i);
        for (int j = 0; j < task.size; ++j) {
            group[j] = task.start + j;
        }
    }
    mesh.update_num();
    return task_list[MPI::rank];
}
