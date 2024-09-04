#include "solver/solver.h"

using namespace MESO;

MESO::MPI::MPI_TaskObject MESO::MPI::DVS_partition(fvmMesh::Mesh &mesh) {
    auto task_list = MPI::get_task_distribution(mesh.NCELL);
    mesh.cell_names.resize(MPI::processor_num);
    mesh.cell_groups.resize(MPI::processor_num);
    mesh.cell_partition_groups.resize(MPI::processor_num);
    for (int i = 0; i < MPI::processor_num; ++i) {
        auto &task = task_list[i];
        auto &group = mesh.cell_groups[i];
        auto &partition = mesh.cell_partition_groups[i];
        group.resize(task.size);
        partition.resize(task.size);
        mesh.cell_names[i] = "node-" + std::to_string(i);
        for (int j = 0; j < task.size; ++j) {
            const int id = task.start + j;
            auto &cell = mesh.cells[id];
            cell.partition_id = i;
            cell.partition_cell_id = j;
            group[j] = id;
            partition[j] = id;
        }
    }
    mesh.update_mesh_params();
    return task_list[MPI::rank];
}
