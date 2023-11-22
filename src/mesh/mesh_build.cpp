#include "mesh.h"

/// node set map
std::unordered_map<int, std::vector<std::vector<int>>> node_set_map = {
        {GEOM::TRIA,  {
                              {0, 1}, // face 1
                              {1, 2}, // face 2
                              {2, 0}, // face 3
                      }},
        {GEOM::QUAD,  {
                              {0, 1}, // face 1
                              {1, 2}, // face 2
                              {2, 3}, // face 3
                              {3, 0}, // face 4
                      }},
        {GEOM::HEXAH, {
                              {0, 1, 2, 3}, // face 1
                              {0, 1, 5, 4}, // face 2
                              {1, 2, 6, 5}, // face 3
                              {2, 3, 7, 6}, // face 4
                              {0, 3, 7, 4}, // face 5
                              {4, 5, 6, 7}, // face 6
                      }},
        {GEOM::PRISM, {
                              {0, 1, 2, 3}, // face 1
                              {0, 1, 4},    // face 2
                              {1, 2, 4},    // face 3
                              {2, 3, 4},    // face 4
                              {3, 0, 4},    // face 5
                      }},
        {GEOM::TETRA, {
                              {0, 1, 2}, // face 1
                              {0, 1, 3}, // face 2
                              {0, 2, 3}, // face 3
                              {1, 2, 3}, // face 4
                      }},
        {GEOM::PYRAM, {
                              {0, 1, 4, 3},   // face 1
                              {1, 2, 5, 4},   // face 2
                              {0, 3, 5, 2},   // face 3
                              {0, 1, 2},      // face 4
                              {3, 4, 5},      // face 5
                      }},
};

/// 内部函数声明
TP_mesh
void build_face_2d(mesh_type &mesh);
/// explicit init
TP_func void build_face_2d<MESH::ListMesh>(MESH::ListMesh &mesh);
/// explicit init
TP_func void build_face_2d<MESH::MapMesh>(MESH::MapMesh &mesh);

TP_mesh
void build_face_3d(mesh_type &mesh);
/// explicit init
TP_func void build_face_3d<MESH::ListMesh>(MESH::ListMesh &mesh);
/// explicit init
TP_func void build_face_3d<MESH::MapMesh>(MESH::MapMesh &mesh);

TP_key_mesh
void generate_face(mesh_type &mesh, MESH::Cell<key_type> &cell,
                   int face_type, int face_on_cell_key,
                   std::unordered_map<std::string, key_type> &face_map);

TP_key_mesh
MESH::Face<key_type> &create_face(mesh_type &mesh, MESH::Cell<key_type> &on_cell, key_vector &node_set,
                                  int face_type, int face_on_cell_key);

/// 函数实现
/// for ListMesh only
void MESH::ListMesh::build_geom() {
    /// cell
    for (auto &cell : CELLS) {
        cell.position = GEOM::cell_position(cell, *this);
        cell.position_square = cell.position * cell.position;
        update_max_discrete_velocity(cell.position.magnitude());
        cell.volume = GEOM::cell_volume(cell, *this);
        update_min_mesh_size(pow(cell.volume, 1.0 / dimension()));
    }
    info_println(" - build geom: cell - ok.");
    if (type == MESH_TYPE_NORMAL) {
        /// face
        for (auto &face : FACES) {
            face.position = GEOM::face_position(face, *this);
            face.area = GEOM::face_area<int, MESH::ListMesh>(face, *this);
            GEOM::face_normal_vector(face, *this);
        }
        info_println(" - build geom: face - ok.");
        /// link near cells by face
        std::unordered_map<int, std::vector<int>> face_map;
        for (auto &cell : CELLS) {
            for (auto &face_key : cell.face_key) {
                face_map[face_key].push_back(cell.key);
            }
        }
        for (auto &it : face_map) {
            auto &face_key = it.first;
            for (auto &cell_key : it.second) {
                auto &cell = get_cell(cell_key);
                for (auto &near_key : it.second) {
                    if (cell_key != near_key) cell.near_cell_key.push_back(near_key);
                }
            }
        }
        /// link second near cells by face
        for (auto &cell : CELLS) {
            for (auto &near_cell_key : cell.near_cell_key) {
                auto &near_cell = get_cell(near_cell_key);
                for (auto &second_near_cell_key : near_cell.near_cell_key) {
                    /// don't append itself
                    if (cell.key == second_near_cell_key) continue;
                    /// don't append cell_key that already exists in cell.near_cell_key
                    if (std::find(cell.near_cell_key.begin(), cell.near_cell_key.end(), second_near_cell_key) != cell.near_cell_key.end()) continue;
                    /// don't append the same cell_key in cell.second_near_cell_key
                    if (std::find(cell.second_near_cell_key.begin(), cell.second_near_cell_key.end(), second_near_cell_key) != cell.second_near_cell_key.end()) continue;
                    /// then we can append!
                    cell.second_near_cell_key.push_back(second_near_cell_key);
                }
            }
        }
    }
    note_println("Build geom - ok.");
}

void MESH::ListMesh::build_face() {
    const int D = dimension();
    if (D == 2) build_face_2d(*this);
    else if (D == 3) build_face_3d(*this);
    note_println("Build interface - ok.");
}

TP_func
void build_face_2d<MESH::ListMesh>(MESH::ListMesh &mesh) {
    /// 2D config - LINE only
    std::unordered_map<std::string, int> line_map;
    for (auto &cell : mesh.CELLS) {
        switch (cell.type) {
            case GEOM::TRIA:
                generate_face(mesh, cell, GEOM::LINE, 0, line_map);
                generate_face(mesh, cell, GEOM::LINE, 1, line_map);
                generate_face(mesh, cell, GEOM::LINE, 2, line_map);
                break;
            case GEOM::QUAD:
                generate_face(mesh, cell, GEOM::LINE, 0, line_map);
                generate_face(mesh, cell, GEOM::LINE, 1, line_map);
                generate_face(mesh, cell, GEOM::LINE, 2, line_map);
                generate_face(mesh, cell, GEOM::LINE, 3, line_map);
                break;
            default:
                std::string s = "build_face_2d() caught invalid geom type.";
                warn_println(s);
                throw std::invalid_argument(s);
        }
    }
    /// link face to mark
    for (int mark_key = 0; mark_key < mesh.MARKS.size(); mark_key++) {
        auto &mark = mesh.get_mark(mark_key);
        for (auto &mark_elem : mark.MARK_ELEM) {
            /// 2D - LINE only
            auto &face_key = line_map.at(mark_elem.key);
            auto &face = mesh.get_face(face_key);
            face.mark_key = mark_key;
            mark_elem.face_key = face_key;
        }
    }
}

TP_func
void build_face_3d<MESH::ListMesh>(MESH::ListMesh &mesh) {
    /// 3D config - quad, tria
    std::unordered_map<std::string, int> quad_map, tria_map;
    for (auto &cell : mesh.CELLS) {
        switch (cell.type) {
            case GEOM::HEXAH:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 1, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 2, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 3, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 4, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 5, quad_map);
                break;
            case GEOM::PRISM:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::TRIA, 1, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 2, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 4, tria_map);
                break;
            case GEOM::TETRA:
                generate_face(mesh, cell, GEOM::TRIA, 0, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 1, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 2, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                break;
            case GEOM::PYRAM:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 1, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 2, quad_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 4, tria_map);
                break;
            default:
                std::string s = "build_face_3d() caught invalid geom type.";
                warn_println(s);
                throw std::invalid_argument(s);
        }
    }
    /// link face to mark
    for (int mark_key = 0; mark_key < mesh.MARKS.size(); mark_key++) {
        auto &mark = mesh.get_mark(mark_key);
        for (auto &mark_elem : mark.MARK_ELEM) {
            switch (mark_elem.type) {
                case GEOM::TRIA: {
                    auto &face_key = tria_map.at(mark_elem.key);
                    auto &face = mesh.get_face(face_key);
                    face.mark_key = mark_key;
                    mark_elem.face_key = face_key;
                    break;
                }
                case GEOM::QUAD: {
                    auto &face_key = quad_map.at(mark_elem.key);
                    auto &face = mesh.get_face(face_key);
                    face.mark_key = mark_key;
                    mark_elem.face_key = face_key;
                    break;
                }
                default:
                    warn_println("Caught mark elem with unsupported geom type.");
                    break;
            }
        }
    }
}

/// for MapMesh only
void MESH::MapMesh::build_geom() {
    /// cell
    for (auto &it : CELLS) {
        auto &cell = it.second;
        cell.position = GEOM::cell_position(cell, *this);
        cell.position_square = cell.position * cell.position;
        update_max_discrete_velocity(cell.position.magnitude());
        cell.volume = GEOM::cell_volume(cell, *this);
        update_min_mesh_size(pow((Pi_over_2 / dimension()) * cell.volume, 1.0 / dimension()));
    }
    info_println(" - build geom: cell - ok.");
    if (type == MESH_TYPE_NORMAL) {
        /// face
        for (auto &it : FACES) {
            auto &face = it.second;
            face.position = GEOM::face_position(face, *this);
            face.area = GEOM::face_area<std::string, MESH::MapMesh>(face, *this);
            GEOM::face_normal_vector(face, *this);
        }
        info_println(" - build geom: face - ok.");
        /// link near cells by face
        std::unordered_map<std::string, string_vector > face_map;
        for (auto &it : CELLS) {
            auto &cell = it.second;
            for (auto &face_key : cell.face_key) {
                face_map[face_key].push_back(cell.key);
            }
        }
        for (auto &it : face_map) {
            auto &face_key = it.first;
            for (auto &cell_key : it.second) {
                auto &cell = get_cell(cell_key);
                for (auto &near_key : it.second) {
                    if (cell_key != near_key) cell.near_cell_key.push_back(near_key);
                }
            }
        }
    }
    note_println("Build geom - ok.");
}

void MESH::MapMesh::build_face() {
    const int D = dimension();
    if (D == 2) build_face_2d(*this);
    else if (D == 3) build_face_3d(*this);
    note_println("Build interface - ok.");
}

TP_func
void build_face_2d<MESH::MapMesh>(MESH::MapMesh &mesh) {
    /// 2D config - LINE only
    std::unordered_map<std::string, std::string> line_map;
    for (auto &it : mesh.CELLS) {
        auto &cell = it.second;
        switch (cell.type) {
            case GEOM::TRIA:
                generate_face(mesh, cell, GEOM::LINE, 0, line_map);
                generate_face(mesh, cell, GEOM::LINE, 1, line_map);
                generate_face(mesh, cell, GEOM::LINE, 2, line_map);
                break;
            case GEOM::QUAD:
                generate_face(mesh, cell, GEOM::LINE, 0, line_map);
                generate_face(mesh, cell, GEOM::LINE, 1, line_map);
                generate_face(mesh, cell, GEOM::LINE, 2, line_map);
                generate_face(mesh, cell, GEOM::LINE, 3, line_map);
                break;
            default:
                std::string s = "build_face_2d() caught invalid geom type.";
                warn_println(s);
                throw std::invalid_argument(s);
        }
    }
    /// link face to mark
    for (auto &it : mesh.MARKS) {
        auto &mark = it.second;
        for (auto &it2 : mark.MARK_ELEM) {
            auto &mark_elem = it2.second;
            /// 2D - LINE only
            auto &face_key = line_map.at(mark_elem.key);
            auto &face = mesh.get_face(face_key);
            face.mark_key = mark.name;
            mark_elem.face_key = face_key;
        }
    }
}

TP_func
void build_face_3d<MESH::MapMesh>(MESH::MapMesh &mesh) {
    /// 3D config - quad, tria
    std::unordered_map<std::string, std::string> quad_map, tria_map;
    for (auto &it : mesh.CELLS) {
        auto &cell = it.second;
        switch (cell.type) {
            case GEOM::HEXAH:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 1, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 2, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 3, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 4, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 5, quad_map);
                break;
            case GEOM::PRISM:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::TRIA, 1, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 2, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 4, tria_map);
                break;
            case GEOM::TETRA:
                generate_face(mesh, cell, GEOM::TRIA, 0, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 1, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 2, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                break;
            case GEOM::PYRAM:
                generate_face(mesh, cell, GEOM::QUAD, 0, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 1, quad_map);
                generate_face(mesh, cell, GEOM::QUAD, 2, quad_map);
                generate_face(mesh, cell, GEOM::TRIA, 3, tria_map);
                generate_face(mesh, cell, GEOM::TRIA, 4, tria_map);
                break;
            default:
                std::string s = "build_face_3d() caught invalid geom type.";
                warn_println(s);
                throw std::invalid_argument(s);
        }
    }
    /// link face to mark
    for (auto &it : mesh.MARKS) {
        auto &mark = it.second;
        for (auto &it2 : mark.MARK_ELEM) {
            auto &mark_elem = it2.second;
            switch (mark_elem.type) {
                case GEOM::TRIA: {
                    auto &face_key = tria_map.at(mark_elem.key);
                    auto &face = mesh.get_face(face_key);
                    face.mark_key = mark.name;
                    mark_elem.face_key = face_key;
                    break;
                }
                case GEOM::QUAD: {
                    auto &face_key = quad_map.at(mark_elem.key);
                    auto &face = mesh.get_face(face_key);
                    face.mark_key = mark.name;
                    mark_elem.face_key = face_key;
                    break;
                }
                default:
                    warn_println("Caught mark elem with unsupported geom type.");
                    break;
            }
        }
    }
}

/// for both Mesh

TP_key_mesh
void generate_face(mesh_type &mesh, MESH::Cell<key_type> &cell,
                   int face_type, int face_on_cell_key, std::unordered_map<std::string, key_type> &face_map) {
    key_vector node_set;
    for (int it : node_set_map[cell.type][face_on_cell_key]) node_set.push_back(cell.node_key[it]);
    node_set.shrink_to_fit();
    /// generate interface key in map.
    std::string key = GEOM::generate_face_key<key_type>(node_set);
    /// check whether exist
    if (face_map.find(key) != face_map.end()) {
        /// exist
        auto &face_key = face_map[key];
        auto &face = mesh.get_face(face_key);
        /// bind face to cell
        face.inv_cell_key = cell.key;
        face.inv_cell_face = face_on_cell_key;
        cell.face_key[face_on_cell_key] = face.key;
        return;
    }
    /// not exist - create face
    auto &face = create_face<key_type, mesh_type>(mesh, cell, node_set, face_type, face_on_cell_key);
    cell.face_key[face_on_cell_key] = face.key;
    /// register face in map
    face_map[key] = face.key;
}

TP_func MESH::Face<int> &create_face(MESH::ListMesh &mesh,
                                     MESH::Cell<int> &on_cell, std::vector<int> &node_set,
                                     int face_type, int face_on_cell_key) {
    int face_key = mesh.face_num();
    mesh.FACES.emplace_back(face_key, face_type, node_set);
    auto &face = mesh.get_face(face_key);
    face.on_cell_key = face.inv_cell_key = on_cell.key;
    face.on_cell_face = face.inv_cell_face = face_on_cell_key;
    face.node_key.resize(node_set.size());
    std::copy(node_set.begin(), node_set.end(), face.node_key.begin());
    return face;
}

TP_func MESH::Face<std::string> &create_face(MESH::MapMesh &mesh,
                                             MESH::Cell<std::string> &on_cell, std::vector<std::string> &node_set,
                                     int face_type, int face_on_cell_key) {
    std::string face_key = std::to_string(mesh.face_num());
    mesh.FaceKey.emplace_back(face_key);
    mesh.FACES.emplace(face_key, MESH::Face<std::string>(face_key, face_type, node_set));
    auto &face = mesh.get_face(face_key);
    face.on_cell_key = face.inv_cell_key = on_cell.key;
    face.on_cell_face = face.inv_cell_face = face_on_cell_key;
    face.node_key.resize(node_set.size());
    std::copy(node_set.begin(), node_set.end(), face.node_key.begin());
    return face;
}
