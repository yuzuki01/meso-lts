#include "mesh.h"

/// su2 mesh
TP_func
void su2Reader::parse(MESH::StaticMesh &mesh) {
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        /// skip empty
        if (data.empty() || std::strncmp(data[0].c_str(), "%", 1) == 0) continue;
        switch (read_case) {
            case 0:
                if (str_vec_cmp(data, ndime_mark)) {
                    mesh.NDIME = stoi(data[1]);
                    break;
                }
                if (str_vec_cmp(data, npoin_mark)) {
                    mesh.NPOIN = stoi(data[1]);
                    read_case = 1;
                    break;
                }
                if (str_vec_cmp(data, nelem_mark)) {
                    mesh.NELEM = stoi(data[1]);
                    read_case = 2;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                break;
            case 1:     // node
                if (str_vec_cmp(data, nelem_mark)) {
                    mesh.NELEM = stoi(data[1]);
                    read_case = 2;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                {
                    Vec3D position{0.0,0.0,0.0};
                    const int len = int(data.size());
                    if (len > 0) position.x = stod(data[0]);
                    if (len > 1) position.y = stod(data[1]);
                    if (len > 2 && mesh.NDIME == 3) position.z = stod(data[2]);
                    mesh.NODES.emplace_back(mesh.NODES.size(), position);
                }
                break;
            case 2:     // cell
                if (str_vec_cmp(data, npoin_mark)) {
                    mesh.NPOIN = stoi(data[1]);
                    read_case = 1;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                {
                    int type = stoi(data[0]);
                    const int node_num = GEOM::node_num(type);
                    std::vector<int> node_set(node_num, -1);
                    for (int k = 0; k < node_num; k++) {
                        node_set[k] = stoi(data[k + 1]);
                    }
                    mesh.CELLS.emplace_back(mesh.CELLS.size(), type, node_set);
                }
                break;
            case 3:     // mark
                if (str_vec_cmp(data, markertag_mark)) {
                    int nelem = 0;
                    std::string mark_name = data[1];
                    i++;
                    while (true) {
                        data = split(lines[i]);
                        if (str_vec_cmp(data, markerelems_mark)) {
                            nelem = stoi(data[1]);
                            break;
                        }
                        if (i >= lines.size()) break;
                        i++;
                    }
                    mesh.MARKS.emplace_back(mark_name, nelem);
                    break;
                }
                {
                    auto &mark = mesh.MARKS.back();
                    mark.MARK_ELEM.emplace_back(data);
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
void su2Reader::parse(MESH::MapMesh &mesh) {
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        /// skip empty
        if (data.empty() || std::strncmp(data[0].c_str(), "%", 1) == 0) continue;
        switch (read_case) {
            case 0:
                if (str_vec_cmp(data, ndime_mark)) {
                    mesh.NDIME = stoi(data[1]);
                    break;
                }
                if (str_vec_cmp(data, npoin_mark)) {
                    mesh.NPOIN = stoi(data[1]);
                    read_case = 1;
                    break;
                }
                if (str_vec_cmp(data, nelem_mark)) {
                    mesh.NELEM = stoi(data[1]);
                    read_case = 2;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                break;
            case 1:     // node
                if (str_vec_cmp(data, nelem_mark)) {
                    mesh.NELEM = stoi(data[1]);
                    read_case = 2;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                {
                    Vec3D position{0.0,0.0,0.0};
                    const int len = int(data.size());
                    if (len > 0) position.x = stod(data[0]);
                    if (len > 1) position.y = stod(data[1]);
                    if (len > 2 && mesh.NDIME == 3) position.z = stod(data[2]);
                    int key = int(mesh.NodeKey.size());
                    mesh.NodeKey.push_back(key);
                    mesh.NODES.emplace(key, MESH::Node(key, position));
                }
                break;
            case 2:     // cell
                if (str_vec_cmp(data, npoin_mark)) {
                    mesh.NPOIN = stoi(data[1]);
                    read_case = 1;
                    break;
                }
                if (str_vec_cmp(data, nmark_mark)) {
                    mesh.NMARK = stoi(data[1]);
                    read_case = 3;
                    break;
                }
                {
                    int type = stoi(data[0]);
                    const int node_num = GEOM::node_num(type);
                    key_vector node_set(node_num, -1);
                    for (int k = 0; k < node_num; k++) {
                        node_set[k] = stoi(data[k + 1]);
                    }
                    int key = int(mesh.CellKey.size());
                    mesh.CellKey.push_back(key);
                    mesh.CELLS.emplace(key, MESH::Cell(key, type, node_set));
                }
                break;
            case 3:     // mark
                debug_println("read " + lines[i]);
                if (str_vec_cmp(data, markertag_mark)) {
                    int nelem = 0;
                    std::string mark_name = data[1];
                    i++;
                    while (true) {
                        data = split(lines[i]);
                        if (str_vec_cmp(data, markerelems_mark)) {
                            nelem = stoi(data[1]);
                            break;
                        }
                        if (i >= lines.size()) break;
                        i++;
                    }
                    int key = int(mesh.MarkKey.size());
                    mesh.MarkKey.push_back(key);
                    mesh.MARKS.emplace(key, MESH::Mark(mark_name, nelem));
                    break;
                }
                {
                    auto &mark_key = mesh.MarkKey.back();
                    auto &mark = mesh.MARKS.at(mark_key);
                    mark.MARK_ELEM.emplace_back(data);
                }
                break;
            default:
                break;
        }
    }
    logger << "Parsed file: " << path;
    logger.info();
}
