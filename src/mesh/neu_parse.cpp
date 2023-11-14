#include "mesh.h"

TP_func
void neuReader::parse(MESH::MapMesh &mesh) {
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        /// skip empty
        if (data.empty()) continue;
        // debug_println("read: " + lines[i]);
        switch (read_case) {
            case 0: // find CONTROL INFO
                if (str_vec_cmp(data, "CONTROL INFO")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (str_vec_cmp(data, "NUMNP NELEM NGRPS NBSETS NDFCD NDFVL")) {
                            data = split(lines[++i]);
                            mesh.set_mesh_params(data);
                        }
                    }
                    read_case = 1;
                }
                break;
            case 1: // find NODAL COORDINATES
                if (str_vec_cmp(data, "NODAL COORDINATES")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        auto tec_id = data[0];
                        mesh.NodeKey.push_back(tec_id);
                        mesh.NODES.emplace(tec_id, MESH::Node<std::string>(data));
                    }
                    read_case = 2;
                }
                break;
            case 2: // find ELEMENTS/CELLS
                if (str_vec_cmp(data, "ELEMENTS/CELLS")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (int(data.size()) < 3 + stoi(data[2])) {
                            string_vector next_line = split(lines[++i]);
                            for (auto &it : next_line) data.push_back(it);
                        }
                        auto tec_id = data[0];
                        mesh.CellKey.push_back(tec_id);
                        mesh.CELLS.emplace(tec_id, MESH::Cell<std::string>(data));
                    }
                    read_case = 3;
                }
                break;
            case 3: // find BOUNDARY CONDITIONS
                if (str_vec_cmp(data, "BOUNDARY CONDITIONS")) {
                    data = split(lines[++i]);
                    if (data.empty()) continue;
                        auto mark_name = data[0];
                        mesh.MarkKey.push_back(mark_name);
                        mesh.MARKS.emplace(mark_name, MESH::Mark<std::string>(data));
                        auto &mark = mesh.get_mark(mark_name);
                        while (true) {
                            data = split(lines[++i]);
                            if (data.empty()) continue;
                            if (data[0] == end_mark) break;
                            auto mark_elem_key = std::to_string(mark.MarkElemKey.size());
                            mark.MarkElemKey.push_back(mark_elem_key);
                            mark.MARK_ELEM.emplace(mark_elem_key, MESH::MarkElem<std::string>(data));
                        }
                }
                break;
            default:
                break;
        }
    }
    logger << "Parsed file: " << path;
    logger.info();
}

TP_func
void neuReader::parse(MESH::ListMesh &mesh) {
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        /// skip empty
        if (data.empty()) continue;
        // debug_println("read: " + lines[i]);
        switch (read_case) {
            case 0: // find CONTROL INFO
                if (str_vec_cmp(data, "CONTROL INFO")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (str_vec_cmp(data, "NUMNP NELEM NGRPS NBSETS NDFCD NDFVL")) {
                            data = split(lines[++i]);
                            mesh.set_mesh_params(data);
                        }
                    }
                    read_case = 1;
                }
                break;
            case 1: // find NODAL COORDINATES
                if (str_vec_cmp(data, "NODAL COORDINATES")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        mesh.NODES.emplace_back(data);
                    }
                    read_case = 2;
                }
                break;
            case 2: // find ELEMENTS/CELLS
                if (str_vec_cmp(data, "ELEMENTS/CELLS")) {
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        if (int(data.size()) < 3 + stoi(data[2])) {
                            string_vector next_line = split(lines[++i]);
                            for (auto &it : next_line) data.push_back(it);
                        }
                        mesh.CELLS.emplace_back(data);
                    }
                    read_case = 3;
                }
                break;
            case 3: // find BOUNDARY CONDITIONS
                if (str_vec_cmp(data, "BOUNDARY CONDITIONS")) {
                    data = split(lines[++i]);
                    if (data.empty()) continue;
                    mesh.MARKS.emplace_back(data);
                    auto &mark = mesh.MARKS.back();
                    while (true) {
                        data = split(lines[++i]);
                        if (data.empty()) continue;
                        if (data[0] == end_mark) break;
                        mark.MARK_ELEM.emplace_back(data);
                    }
                }
                break;
            default:
                break;
        }
    }
    logger << "Parsed file: " << path;
    logger.info();
}
