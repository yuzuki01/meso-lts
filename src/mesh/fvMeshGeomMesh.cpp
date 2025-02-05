#include "mesh/mesh.h"

namespaceMesoMesh


GeomMesh::GeomMesh(const FileReader::BasicReader &reader)
        : name(reader.name()) {
    parseGambit(reader.read_lines());
}

void GeomMesh::parseGambit(const List<String> &lines) {
    int line_size = int(lines.size());
    StringList data;
    for (int i = 0; i < line_size; ++i) {
        data = MESO::Utils::split(lines[i]);
        if (data[0] == "NUMNP") {
            data = MESO::Utils::split(lines[++i]);
            this->NNODE = std::stoi(data[0]);
            this->NCELL = std::stoi(data[1]);
            this->NZONE = std::stoi(data[2]);
            zones_.reserve(NZONE);
            this->NMARK = std::stoi(data[3]);
            this->NDFCD = std::stoi(data[4]);
            this->NDFVL = std::stoi(data[5]);
            continue;
        }
        if (data[0] == "NODAL" and data[1] == "COORDINATES") {
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse node
                int node_id = std::stoi(data[0]) - 1;
                double nx = 0.0, ny = 0.0, nz = 0.0;
                if (data.size() >= 3) {
                    nx = std::stod(data[1]);
                    ny = std::stod(data[2]);
                }
                if (data.size() >= 4) {
                    nz = std::stod(data[3]);
                }
                nodes_.emplace_back(*this,
                                    node_id,
                                    Vector(nx, ny, nz));
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings nodes, NDODE=" << this->NNODE << std::endl;
            continue;
        }
        if (data[0] == "ELEMENTS/CELLS") {
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse cell
                int cell_id = std::stoi(data[0]) - 1;
                int geom_type = std::stoi(data[1]);
                int node_count = std::stoi(data[2]);
                int data_len = int(data.size());
                List<ObjectId> node_list;
                for (int j = 3; j < data_len; ++j) {
                    int node_id = std::stoi(data[j]) - 1;
                    node_list.push_back(node_id);
                }
                if (node_list.size() < node_count) {
                    ++i;
                    data = MESO::Utils::split(lines[i]);
                    for (const auto &it: data) {
                        int node_id = std::stoi(it) - 1;
                        node_list.push_back(node_id);
                    }
                }
                cells_.emplace_back(*this, cell_id, geom_type, node_list);
                ++i;
                /// generate face

            }
            logger.debug << "Read " << name << " for GeomMesh, NCELL=" << this->NCELL << std::endl;
            continue;
        }
        if (data[0] == "ELEMENT" and data[1] == "GROUP") {
            ++i;
            data = MESO::Utils::split(lines[i]);
            int group_id = std::stoi(data[1]) - 1;
            // int group_nelem = std::stoi(data[3]);
            ++i;
            data = MESO::Utils::split(lines[i]);
            std::string group_name = data[0];
            zones_.emplace_back(*this, group_name);
            i += 2;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse
                for (const auto &it: data) {
                    int cell_id = std::stoi(it) - 1;
                    auto &cell = cells_[cell_id];
                    zones_[group_id].group().push_back(cell_id);
                }
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings group: " << group_name << std::endl;
            continue;
        }
        if (data[0] == "BOUNDARY" and data[1] == "CONDITIONS") {
            ++i;
            data = MESO::Utils::split(lines[i]);
            std::string bc_name = data[0];
            int face_num = stoi(data[2]);
            int group_id = int(marks_.size());
            marks_.emplace_back(*this, bc_name);
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse bset
                int cell_id = std::stoi(data[0]) - 1;
                int face_on_cell_id = std::stoi(data[2]) - 1;
                auto &cell = cells_[cell_id];
                auto &face = faces_[cell.face_id[face_on_cell_id]];
                face.group_id = group_id;
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings mark: " << bc_name << std::endl;
            continue;
        }
    }
    for (auto &face: this->faces) {
        this->face_groups[face.group_id].push_back(face.id);
    }
}


void GeomMesh::parseGambitHead(Label &i, const Label &size,
                               const StringList &lines) {
    StringList data;
    for (; i < size; i++) {
        data = Utils::split(lines[i]);
        if (data[0] == "NUMNP") {
            data = Utils::split(lines[++i]);
            NNODE = std::stoi(data[0]);
            NCELL = std::stoi(data[1]);
            NZONE = std::stoi(data[2]);
            zones_.reserve(NZONE);
            NMARK = std::stoi(data[3]);
            NDFCD = std::stoi(data[4]);
            NDFVL = std::stoi(data[5]);
            break;
        }
    }
    logger.debug << "Read mesh head - Done" << std::endl;
}

void GeomMesh::parseGambitPoint(Label &i, const Label &size,
                                const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        if (data[0] == "NODAL" and data[1] == "COORDINATES") {
            ++i;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                // node
                Label ni = std::stoi(data[0]) - 1;
                Scalar nx = 0, ny = 0, nz = 0;
                if (data.size() >= 3) {
                    nx = std::stod(data[1]);
                    ny = std::stod(data[2]);
                }
                if (data.size() >= 4) {
                    nz = std::stod(data[3]);
                }
                nodes_.emplace_back(*this, ni,
                                    Vector(nx, ny, nz));
                ++i;
            }
        }
        if (data[0] == "ENDOFSECTION") break;
    }
    logger.debug << "Read mesh points - Done" << std::endl;
}

void GeomMesh::parseGambitCell(Label &i, const Label &size,
                               const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        if (data[0] == "ELEMENTS/CELLS") {
            ++i;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                // cell
                ObjectId ci = std::stoi(data[0]) - 1;
                ObjectType geomType = std::stoi(data[1]);
                Label nNum = std::stoi(data[2]);
                auto dataLen = static_cast<Label>(data.size());
                List<ObjectId> node_list;
                for (Label j = 3; j < dataLen; ++j) {
                    ObjectId ni = std::stoi(data[j]) - 1;
                    node_list.push_back(ni);
                }
                if (node_list.size() < nNum) {
                    ++i;
                    data = Utils::split(lines[i]);
                    for (const auto& it : data) {
                        ObjectId ni = std::stoi(it) - 1;
                        node_list.push_back(ni);
                    }
                }
                cells_.emplace_back(*this, ci,
                                    geomType, node_list);
                ++i;
            }
            if (data[0] == "ENDOFSECTION") break;
        }
    }
    logger.debug << "Read mesh cells - Done" << std::endl;
}


/**
 *  ================================================
 *  ------------------ Interfaces ------------------
 *  ================================================
 **/

const Node &GeomMesh::node(const ObjectId &id) const {
    return nodes_[id];
}

const List<Node> &GeomMesh::nodes() const {
    return nodes_;
}

const Cell &GeomMesh::cell(const ObjectId &id) const {
    return cells_[id];
}

const List<Cell> &GeomMesh::cells() const {
    return cells_;
}
