#include "mesh.h"

TP_func MeshWriter<MESH::StaticMesh>::MeshWriter(std::string _path, MESH::StaticMesh &_mesh)
        : path(std::move(_path)), mesh(_mesh) {
    fp.open(path, std::ios::out | std::ios::trunc);
    is_file_open = fp.is_open();
    if (!is_file_open) {
        warn_println("Cannot write to file: " + path);
        fp.close();
        return;
    }
}

TP_func MeshWriter<MESH::MapMesh>::MeshWriter(std::string _path, MESH::MapMesh &_mesh)
        : path(std::move(_path)), mesh(_mesh) {
    fp.open(path, std::ios::out | std::ios::trunc);
    is_file_open = fp.is_open();
    if (!is_file_open) {
        warn_println("Cannot write to file: " + path);
        fp.close();
        return;
    }
}

TP_func bool MeshWriter<MESH::StaticMesh>::is_open() const {
    return is_file_open;
}

TP_func bool MeshWriter<MESH::MapMesh>::is_open() const {
    return is_file_open;
}

TP_func void MeshWriter<MESH::StaticMesh>::close() {
    fp.close();
}

TP_func void MeshWriter<MESH::MapMesh>::close() {
    fp.close();
}

TP_func void MeshWriter<MESH::StaticMesh>::write_head(const string_vector &values) {
    string_vector value_names;
    const int D = mesh.dimension();
    if (D == 2) value_names = {"X", "Y"};
    if (D == 3) value_names = {"X", "Y", "Z"};
    for (auto &it : values) value_names.push_back(it);
    const int var_num = int(value_names.size());
    /// write head
    fp << "VARIABLES = ";
    for (auto &it : value_names) fp << "\"" << it << "\",";
    fp << std::endl;
    fp << "ZONE N=" << mesh.node_num() << ", E=" << mesh.cell_num() << ", VARLOCATION=([1-" << D
       << "]=NODAL, ["
       << ((var_num == 1 + D) ?
           std::to_string(D + 1) :
           (std::to_string(D + 1) + "-" + std::to_string(var_num)))
       << "]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE=" << (D == 2 ? "FEQUADRILATERAL" : "FEBRICK")
       << std::endl;
}

TP_func void MeshWriter<MESH::MapMesh>::write_head(const string_vector &values) {
    string_vector value_names;
    const int D = mesh.dimension();
    if (D == 2) value_names = {"X", "Y"};
    if (D == 3) value_names = {"X", "Y", "Z"};
    for (auto &it : values) value_names.push_back(it);
    const int var_num = int(value_names.size());
    /// write head
    fp << "VARIABLES = ";
    for (auto &it : value_names) fp << "\"" << it << "\",";
    fp << std::endl;
    fp << "ZONE N=" << mesh.node_num() << ", E=" << mesh.cell_num() << ", VARLOCATION=([1-" << D
       << "]=NODAL, ["
       << ((var_num == 1 + D) ?
           std::to_string(D + 1) :
           (std::to_string(D + 1) + "-" + std::to_string(var_num)))
       << "]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE=" << (D == 2 ? "FEQUADRILATERAL" : "FEBRICK")
       << std::endl;
}

TP_func void MeshWriter<MESH::StaticMesh>::write_node() {
    int count;
    fp << std::endl;
    // write node
    count = 0;
    for (auto &node : mesh.NODES) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }

    count = 0;
    for (auto &node : mesh.NODES) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }

    if (mesh.dimension() == 3) {
        count = 0;
        for (auto &node : mesh.NODES) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }

    fp << std::endl;
}

TP_func void MeshWriter<MESH::MapMesh>::write_node() {
    int count;
    fp << std::endl;
    // write node
    count = 0;
    for (auto &it : mesh.NODES) {
        auto &node = it.second;
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }

    count = 0;
    for (auto &it : mesh.NODES) {
        auto &node = it.second;
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }

    if (mesh.dimension() == 3) {
        count = 0;
        for (auto &it : mesh.NODES) {
            auto &node = it.second;
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }

    fp << std::endl;
}

TP_func void MeshWriter<MESH::StaticMesh>::write_data(const std::vector<double> &data) {
    int count = 0;
    fp << std::endl << "## cell value" << std::endl;
    for (auto &cell : mesh.CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << data[cell.key];
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
}

TP_func void MeshWriter<MESH::StaticMesh>::write_geom() {
    fp << std::endl << "# geom" << std::endl;
    for (auto &cell : mesh.CELLS) {
        switch (cell.type) {
            case GEOM::TRIA:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1;
                break;
            case GEOM::QUAD:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1;
                break;
            case GEOM::TETRA:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1;
                break;
            case GEOM::HEXAH:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[5]].key + 1
                   << " " << mesh.NODES[cell.node_key[6]].key + 1 << " " << mesh.NODES[cell.node_key[7]].key + 1;
                break;
            case GEOM::PRISM:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1;
                break;
            case GEOM::PYRAM:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1
                   << " " << mesh.NODES[cell.node_key[5]].key + 1 << " " << mesh.NODES[cell.node_key[5]].key + 1;
                break;
            default:
                fp << " unsupported type.";
        }
        fp << std::endl;
    }
}

TP_func void MeshWriter<MESH::MapMesh>::write_geom() {
    std::unordered_map<int, int> node_id, cell_id;
    for (int i = 0; i < mesh.node_num(); i++) node_id[mesh.NodeKey[i]] = i + 1;
    fp << std::endl << "# geom" << std::endl;
    for (auto &it : mesh.CELLS) {
        auto &cell_key = it.first;
        auto &cell = it.second;
        switch (cell.type) {
            case GEOM::TRIA:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[2]] + 1;
                break;
            case GEOM::QUAD:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[3]] + 1;
                break;
            case GEOM::TETRA:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[2]] + 1
                   << " " << node_id[cell.node_key[3]] + 1 << " " << node_id[cell.node_key[3]] + 1
                   << " " << node_id[cell.node_key[3]] + 1 << " " << node_id[cell.node_key[3]] + 1;
                break;
            case GEOM::HEXAH:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[3]] + 1
                   << " " << node_id[cell.node_key[4]] + 1 << " " << node_id[cell.node_key[5]] + 1
                   << " " << node_id[cell.node_key[6]] + 1 << " " << node_id[cell.node_key[7]] + 1;
                break;
            case GEOM::PRISM:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[3]] + 1
                   << " " << node_id[cell.node_key[4]] + 1 << " " << node_id[cell.node_key[4]] + 1
                   << " " << node_id[cell.node_key[4]] + 1 << " " << node_id[cell.node_key[4]] + 1;
                break;
            case GEOM::PYRAM:
                fp << " " << node_id[cell.node_key[0]] + 1 << " " << node_id[cell.node_key[1]] + 1
                   << " " << node_id[cell.node_key[2]] + 1 << " " << node_id[cell.node_key[2]] + 1
                   << " " << node_id[cell.node_key[3]] + 1 << " " << node_id[cell.node_key[4]] + 1
                   << " " << node_id[cell.node_key[5]] + 1 << " " << node_id[cell.node_key[5]] + 1;
                break;
            default:
                fp << " unsupported type.";
        }
        fp << std::endl;
    }
}
