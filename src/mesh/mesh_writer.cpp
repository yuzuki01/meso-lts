#include "mesh.h"

TP_func MeshWriter<int>::MeshWriter(std::string _path, MESH::Mesh<int> &_mesh)
        : path(std::move(_path)), mesh(_mesh) {
    fp.open(path, std::ios::out | std::ios::trunc);
    is_file_open = fp.is_open();
    if (!is_file_open) {
        warn_println("Cannot write to file: " + path);
        fp.close();
        return;
    }
}

TP_func MeshWriter<std::string>::MeshWriter(std::string _path, MESH::Mesh<std::string> &_mesh)
        : path(std::move(_path)), mesh(_mesh) {
    fp.open(path, std::ios::out | std::ios::trunc);
    is_file_open = fp.is_open();
    if (!is_file_open) {
        warn_println("Cannot write to file: " + path);
        fp.close();
        return;
    }
}

TP_func bool MeshWriter<int>::is_open() const {
    return is_file_open;
}

TP_func bool MeshWriter<std::string>::is_open() const {
    return is_file_open;
}

TP_func void MeshWriter<int>::close() {
    fp.close();
}

TP_func void MeshWriter<std::string>::close() {
    fp.close();
}

TP_func void MeshWriter<int>::write_head(const string_vector &values) {
    string_vector value_names;
    const int D = mesh.dimension();
    if (D == 2) value_names = {"X", "Y"};
    if (D == 3) value_names = {"X", "Y", "Z"};
    for (auto &it : values) value_names.push_back(it);
    const int var_num = value_names.size();
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

TP_func void MeshWriter<std::string>::write_head(const string_vector &values) {
    string_vector value_names;
    const int D = mesh.dimension();
    if (D == 2) value_names = {"X", "Y"};
    if (D == 3) value_names = {"X", "Y", "Z"};
    for (auto &it : values) value_names.push_back(it);
    const int var_num = value_names.size();
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

TP_func void MeshWriter<int>::write_node() {
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

TP_func void MeshWriter<std::string>::write_node() {
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

TP_func void MeshWriter<int>::write_data(const std::vector<double> &data) {
    int count = 0;
    fp << std::endl << "## cell value" << std::endl;
    for (auto & cell : mesh.CELLS) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << data[cell.key];
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
}

TP_func void MeshWriter<std::string>::write_data(const std::unordered_map<std::string, double> &data) {
    int count = 0;
    fp << std::endl << "## cell value" << std::endl;
    for (auto & cell_key : mesh.CellKey) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << data.at(cell_key);
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
}

TP_func void MeshWriter<int>::write_geom() {
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
            case GEOM::BRICK:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[5]].key + 1
                   << " " << mesh.NODES[cell.node_key[7]].key + 1 << " " << mesh.NODES[cell.node_key[6]].key + 1;
                break;
            case GEOM::WEDGE:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1
                   << " " << mesh.NODES[cell.node_key[5]].key + 1 << " " << mesh.NODES[cell.node_key[5]].key + 1;
                break;
            case GEOM::TETRA:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[2]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[3]].key + 1;
                break;
            case GEOM::PYRAM:
                fp << " " << mesh.NODES[cell.node_key[0]].key + 1 << " " << mesh.NODES[cell.node_key[1]].key + 1
                   << " " << mesh.NODES[cell.node_key[3]].key + 1 << " " << mesh.NODES[cell.node_key[2]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1
                   << " " << mesh.NODES[cell.node_key[4]].key + 1 << " " << mesh.NODES[cell.node_key[4]].key + 1;
                break;
            default:
                fp << " unsupported type.";
        }
        fp << std::endl;
    }
}

TP_func void MeshWriter<std::string>::write_geom() {
    std::unordered_map<std::string, int> node_id, cell_id;
    for (int i = 0; i < mesh.node_num(); i++) node_id[mesh.NodeKey[i]] = i + 1;
    fp << std::endl << "# geom" << std::endl;
    for (auto &it : mesh.CELLS) {
        auto &cell_key = it.first;
        auto &cell = it.second;
        switch (cell.type) {
            case GEOM::TRIA:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[2]] << " " << node_id[cell.node_key[2]];
                break;
            case GEOM::QUAD:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[2]] << " " << node_id[cell.node_key[3]];
                break;
            case GEOM::BRICK:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[3]] << " " << node_id[cell.node_key[2]]
                   << " " << node_id[cell.node_key[4]] << " " << node_id[cell.node_key[5]]
                   << " " << node_id[cell.node_key[7]] << " " << node_id[cell.node_key[6]];
                break;
            case GEOM::WEDGE:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[2]] << " " << node_id[cell.node_key[2]]
                   << " " << node_id[cell.node_key[3]] << " " << node_id[cell.node_key[4]]
                   << " " << node_id[cell.node_key[5]] << " " << node_id[cell.node_key[5]];
                break;
            case GEOM::TETRA:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[2]] << " " << node_id[cell.node_key[2]]
                   << " " << node_id[cell.node_key[3]] << " " << node_id[cell.node_key[3]]
                   << " " << node_id[cell.node_key[3]] << " " << node_id[cell.node_key[3]];
                break;
            case GEOM::PYRAM:
                fp << " " << node_id[cell.node_key[0]] << " " << node_id[cell.node_key[1]]
                   << " " << node_id[cell.node_key[3]] << " " << node_id[cell.node_key[2]]
                   << " " << node_id[cell.node_key[4]] << " " << node_id[cell.node_key[4]]
                   << " " << node_id[cell.node_key[4]] << " " << node_id[cell.node_key[4]];
                break;
            default:
                fp << " unsupported type.";
        }
        fp << std::endl;
    }
}
