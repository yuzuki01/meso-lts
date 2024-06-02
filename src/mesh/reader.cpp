#include <utility>

#include "mesh/mesh.h"
#include "mesh/reader.h"


using namespace MESO::Mesh::Reader;


MESO::StringList BasicReader::read_lines() const {
    std::ifstream fp(file);
    if (!fp.is_open()) {
        logger.warn << "Cannot open the file: " << file << std::endl;
        fp.close();
        throw std::invalid_argument("file open error");
    }
    String line;
    StringList lines;
    while (std::getline(fp, line)) {
        /// skip empty line
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    fp.close();
    return lines;
}

MESO::Mesh::Zone Gambit::read() {
    StringList lines = read_lines();
    int line_size = int(lines.size());
    StringList data;
    Zone zone;
    for (int i = 0; i < line_size; ++i) {
        data = MESO::Utils::split(lines[i]);
        if (data[0] == "NNODE") {
            data = MESO::Utils::split(lines[++i]);
            zone.NNODE = std::stoi(data[0]);
            zone.nodes.reserve(zone.NNODE);
            zone.NCELL = std::stoi(data[1]);
            zone.cells.reserve(zone.NCELL);
            zone.NZONE = std::stoi(data[2]);
            zone.cell_names.reserve(zone.NZONE);
            zone.cell_groups.resize(zone.NZONE);
            zone.NMARK = std::stoi(data[3]);
            zone.face_names.reserve(zone.NMARK + 1);
            zone.face_groups.resize(zone.NMARK + 1);
            zone.NDFCD = std::stoi(data[4]);
            zone.NDFVL = std::stoi(data[5]);
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
                zone.nodes.emplace_back(node_id, Vector{nx, ny, nz});
                ++i;
            }
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
                ObjectIdList node_list;
                for (int j = 3; j < data_len; ++j) {
                    node_list.push_back(std::stoi(data[j]) - 1);
                }
                if (node_list.size() < node_count) {
                    ++i;
                    data = MESO::Utils::split(lines[i]);
                    for (const auto& it : data) {
                        node_list.push_back(std::stoi(it));
                    }
                }
                zone.cells.emplace_back(cell_id, geom_type, node_list);
                ++i;
            }
            /// generate face
            zone.generate_face();
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
            zone.cell_names.push_back(group_name);
            i += 2;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse
                for (const auto &it : data) {
                    int cell_id = std::stoi(it) - 1;
                    auto &cell = zone.cells[cell_id];
                    cell.group_id = group_id;
                    zone.cell_groups[group_id].push_back(cell_id);
                }
                ++i;
            }
            continue;
        }
        if (data[0] == "BOUNDARY" and data[1] == "CONDITIONS") {
            ++i;
            data = MESO::Utils::split(lines[i]);
            std::string name = data[0];
            int face_num = stoi(data[2]);
            int group_id = int(zone.face_names.size());
            zone.face_names.push_back(name);
            zone.face_groups[group_id].reserve(face_num);
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse bset
                int cell_id = std::stoi(data[0]) - 1;
                int face_on_cell_id = std::stoi(data[2]) - 1;
                auto &cell = zone.cells[cell_id];
                auto &face = zone.faces[cell.face_id[face_on_cell_id]];
                face.group_id = group_id;
                ++i;
            }
            continue;
        }
    }
    for (auto &face : zone.faces) {
        zone.face_groups[face.group_id].push_back(face.id);
    }
    return zone;
}
