#include "mesh/mesh.h"
#include "mesh/reader.h"


using namespace MESO::fvmMesh::Reader;


MESO::fvmMesh::Mesh Gambit::parse_mesh_strings() {
    StringList lines = read_lines();
    int line_size = int(lines.size());
    StringList data;
    Mesh mesh;
    for (int i = 0; i < line_size; ++i) {
        data = MESO::Utils::split(lines[i]);
        if (data[0] == "NUMNP") {
            data = MESO::Utils::split(lines[++i]);
            mesh.NNODE = std::stoi(data[0]);
            mesh.nodes.reserve(mesh.NNODE);
            mesh.NCELL = std::stoi(data[1]);
            mesh.cells.reserve(mesh.NCELL);
            mesh.NZONE = std::stoi(data[2]);
            mesh.cell_names.reserve(mesh.NZONE);
            mesh.cell_groups.resize(mesh.NZONE);
            mesh.NMARK = std::stoi(data[3]);
            mesh.face_names.reserve(mesh.NMARK + 1);
            mesh.face_groups.resize(mesh.NMARK + 1);
            mesh.NDFCD = std::stoi(data[4]);
            mesh.NDFVL = std::stoi(data[5]);
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
                mesh.nodes.emplace_back(node_id, Vector(nx, ny, nz));
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings nodes, NDODE=" << mesh.NNODE << std::endl;
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
                    int node_id = std::stoi(data[j]) - 1;
                    node_list.push_back(node_id);
                    mesh.nodes[node_id].neighbors.push_back(cell_id);
                }
                if (node_list.size() < node_count) {
                    ++i;
                    data = MESO::Utils::split(lines[i]);
                    for (const auto &it: data) {
                        int node_id = std::stoi(it) - 1;
                        node_list.push_back(node_id);
                        mesh.nodes[node_id].neighbors.push_back(cell_id);
                    }
                }
                mesh.cells.emplace_back(cell_id, geom_type, node_list);
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings cells, NCELL=" << mesh.NCELL << std::endl;
            /// generate face
            mesh.generate_face();
            logger.debug << "mesh - build faces: " << mesh.NFACE << std::endl;
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
            mesh.cell_names.push_back(group_name);
            i += 2;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse
                for (const auto &it: data) {
                    int cell_id = std::stoi(it) - 1;
                    auto &cell = mesh.cells[cell_id];
                    cell.group_id = group_id;
                    mesh.cell_groups[group_id].push_back(cell_id);
                }
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings group: " << group_name << std::endl;
            continue;
        }
        if (data[0] == "BOUNDARY" and data[1] == "CONDITIONS") {
            ++i;
            data = MESO::Utils::split(lines[i]);
            std::string name = data[0];
            int face_num = stoi(data[2]);
            int group_id = int(mesh.face_names.size());
            mesh.face_names.push_back(name);
            mesh.face_groups[group_id].reserve(face_num);
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse bset
                int cell_id = std::stoi(data[0]) - 1;
                int face_on_cell_id = std::stoi(data[2]) - 1;
                auto &cell = mesh.cells[cell_id];
                auto &face = mesh.faces[cell.face_id[face_on_cell_id]];
                face.group_id = group_id;
                ++i;
            }
            logger.debug << "mesh - parse_mesh_strings mark: " << name << std::endl;
            continue;
        }
    }
    for (auto &face: mesh.faces) {
        mesh.face_groups[face.group_id].push_back(face.id);
    }
    return mesh;
}
