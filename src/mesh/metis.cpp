#include "mesh/mesh.h"

using namespace MESO::fvmMesh;

#ifdef _METIS_H_
void Mesh::partition() {
    auto numCells = static_cast<idx_t>(cells.size());
    auto numParts = static_cast<idx_t>(MPI::processor_num);
    if (MPI::processor_num == 1) {
        /// if process_num == 1   mesh has only 1 part.
        cell_partition_groups.resize(numParts, {});
        for (auto &cell : cells) {
            cell.partition_id = MPI::rank;
            cell.partition_cell_id = cell.id;
            cell_partition_groups[0].push_back(cell.id);
        }
        cell_partition_groups[0].shrink_to_fit();

        logger.info << "Partition mesh by METIS into " << numParts << " part(s)." << std::endl;
        return;
    }

    cell_partition_groups.resize(numParts, {});

    typedef std::vector<idx_t> MetisIdList;

    MetisIdList parts(numCells);

    if (MPI::rank == MPI::main_rank) {
        MetisIdList xadj(numCells + 1);
        MetisIdList adjncy;

        xadj[0] = 0;
        for (int i = 0; i < numCells; ++i) {
            Cell& cell = cells[i];
            xadj[i + 1] = xadj[i] + static_cast<int>(cell.neighbors.size());
            for (auto neighbor : cell.neighbors) {
                adjncy.push_back(static_cast<idx_t>(neighbor));
            }
        }

        idx_t objval;
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        METIS_PartGraphKway(&numCells, &numParts, xadj.data(), adjncy.data(),
                            nullptr, nullptr, nullptr, &numParts, nullptr,
                            nullptr, options, &objval, parts.data());
    }
    MPI::Bcast(parts);
    MPI_Barrier(MPI_COMM_WORLD);

    for (idx_t i = 0; i < numCells; ++i) {
        auto &cell = cells[i];
        auto &group = cell_partition_groups[cell.partition_id];
        cell.partition_id = static_cast<int>(parts[i]);
        cell.partition_cell_id = int(group.size());
        group.push_back(cell.id);
    }
    for (auto &group : cell_partition_groups) {
        group.shrink_to_fit();
    }

    logger.info << "Partition mesh by METIS into " << numParts << " part(s)." << std::endl;
}
#endif
