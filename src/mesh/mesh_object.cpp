#include "mesh.h"

using namespace MESH;

/// global vars
std::unordered_map<std::string, int> MESH::MarkTypeID{
        {"interface",       MESH_BC_INTERFACE},
        {"inlet",           MESH_BC_INLET},
        {"outlet",          MESH_BC_OUTLET},
        {"isothermal_wall", MESH_BC_ISOTHERMAL_WALL},
        {"adiabat_wall",    MESH_BC_ABDIABAT_WALL},
        {"symmetry",        MESH_BC_SYMMETRY},
        {"periodic",        MESH_BC_PERIODIC},
};

/// Node
TP_func Node<int>::Node(const int &node_key, const Vec3D &node_position) : key(node_key), position(node_position) {}

TP_func Node<std::string>::Node(const std::string &node_key, const Vec3D &node_position) : key(node_key), position(node_position) {}

/// Cell
TP_func Cell<int>::Cell(const int &cell_key, const int cell_type, const std::vector<int> &node_set) {
    key = cell_key;
    type = cell_type;
    const int node_num = GEOM::node_num(type);
    node_key.resize(node_num);
    std::copy(node_set.begin(), node_set.end(), node_key.begin());
    volume = -1.0;
    face_num = GEOM::face_num(type);
    face_key.resize(face_num, -1);
    position_square = -1.0;
}

TP_func Cell<std::string>::Cell(const std::string &cell_key, const int cell_type, const std::vector<std::string> &node_set) {
    key = cell_key;
    type = cell_type;
    const int node_num = GEOM::node_num(type);
    node_key.resize(node_num);
    std::copy(node_set.begin(), node_set.end(), node_key.begin());
    volume = -1.0;
    face_num = GEOM::face_num(type);
    face_key.resize(face_num, STRING_NULL);
    position_square = -1.0;
}

TP_func Cell<int>::Cell(const Vec3D &particle_velocity, double weight,
                        const int &cell_key, const int &inv_key,
                        int on_layer_cell_id, int on_layer_cell_num) {
    /// structural DVS
    key = cell_key;
    type = -1;
    volume = weight;
    inv_cell_key = inv_key;
    position = particle_velocity;
    on_layer_id = on_layer_cell_id;
    layer_cell_num = on_layer_cell_num;
    face_num = 0;
    position_square = particle_velocity * particle_velocity;
}

TP_func Cell<std::string>::Cell(const Vec3D &particle_velocity, double weight,
                                const std::string &cell_key, const std::string &inv_key,
                                int on_layer_cell_id, int on_layer_cell_num) {
    /// structural DVS
    key = cell_key;
    type = -1;
    volume = weight;
    inv_cell_key = inv_key;
    position = particle_velocity;
    on_layer_id = on_layer_cell_id;
    layer_cell_num = on_layer_cell_num;
    face_num = 0;
    position_square = particle_velocity * particle_velocity;
}

TP_key void Cell<key_type>::shrink_to_fit() {
    node_key.shrink_to_fit();
    face_key.shrink_to_fit();
    near_cell_key.shrink_to_fit();
    second_near_cell_key.shrink_to_fit();
}

/// Face
TP_func Face<int>::Face(const int &face_key, int elem_type, const std::vector<int> &node_key_vec) {
    key = face_key;
    type = elem_type;
    node_key.resize(node_key_vec.size());
    std::copy(node_key_vec.begin(), node_key_vec.end(), node_key.begin());
}

TP_func Face<std::string>::Face(const std::string &face_key, int elem_type, const string_vector &node_key_vec) {
    key = face_key;
    type = elem_type;
    node_key.resize(node_key_vec.size());
    std::copy(node_key_vec.begin(), node_key_vec.end(), node_key.begin());
}

TP_key void Face<key_type>::shrink_to_fit() {
    node_key.shrink_to_fit();
}

TP_func std::string Face<int>::info() const {
    std::stringstream ss;
    ss << "face<key=" << key << ">" << "position=" << position.info();
    return ss.str();
}

TP_func std::string Face<std::string>::info() const {
    std::stringstream ss;
    ss << "face<key=" << key << ">" << "position=" << position.info();
    return ss.str();
}

/// MarkElem
TP_func MarkElem<int>::MarkElem(const std::vector<std::string> &init_strings) {
    type = stoi(init_strings[0]);
    const int node_num = GEOM::node_num(type);
    node_set.resize(node_num, -1);
    for (int i = 0; i < node_num; i++) {
        node_set[i] = stoi(init_strings[i + 1]);
    }
    key = GEOM::generate_face_key(node_set);
}

TP_func MarkElem<std::string>::MarkElem(const std::vector<std::string> &init_strings) {
    type = stoi(init_strings[0]);
    const int node_num = GEOM::node_num(type);
    node_set.resize(node_num, STRING_NULL);
    for (int i = 0; i < node_num; i++) {
        node_set[i] = init_strings[i + 1];
    }
    key = GEOM::generate_face_key(node_set);
}

/// Mark
Mark<int>::Mark(const std::string &mark_name, const int nelem) {
    name = mark_name;
    elem_num = nelem;
    MARK_ELEM.reserve(elem_num);
}

Mark<std::string>::Mark(const std::string &mark_name, const int nelem) {
    name = mark_name;
    elem_num = nelem;
    MARK_ELEM.reserve(elem_num);
}

void Mark<int>::shrink_to_fit() {
    MARK_ELEM.shrink_to_fit();
}

void Mark<std::string>::shrink_to_fit() {
    MarkElemKey.shrink_to_fit();
    MARK_ELEM.rehash(REHASH_SIZE);
}

MarkElem<int> &Mark<int>::get_elem(const int &_key) {
    return MARK_ELEM[_key];
}

MarkElem<std::string> &Mark<std::string>::get_elem(const std::string &_key) {
    return MARK_ELEM.at(_key);
}

void Mark<int>::info() const {
    std::stringstream ss;
    ss << "  Mark<" << name << ">\n"
       << "     - elements: " << elem_num << "\n"
       << "     - type: " << type << "\n"
       << "     - density: " << density << " temperature: " << temperature << "\n"
       << "     - velocity: " << velocity.info();
    info_println(ss.str());
}

void Mark<std::string>::info() const {
    std::stringstream ss;
    ss << "  Mark<" << name << ">\n"
       << "     - elements: " << elem_num << "\n"
       << "     - type: " << type << "\n"
       << "     - density: " << density << " temperature: " << temperature << "\n"
       << "     - velocity: " << velocity.info();
    info_println(ss.str());
}

/// Mesh
///     Basic mesh:
int BasicMesh::dimension() const {
    return NDIME;
}

void BasicMesh::info() {
    if (name == MESH_KEY_NULL) return;
    std::string type_name;
    switch (type) {
        case MESH_TYPE_NORMAL:
            type_name = "normal";
            break;
        case MESH_TYPE_NO_FACE:
            type_name = "no-face";
            break;
        default:
            type_name = "unknown";
    }
    std::stringstream ss;
    note_println("Mesh info:");
    ss << "  file: " << name << "\n"
       << "  type: " << type_name << "\n"
       << "  mesh data:";
    info_println(ss.str());
    data_int_println({"node", "cell", "face", "mark"}, {NPOIN, NELEM, -1, NMARK});
    data_double_println({"max-magnitude", "min-size"}, {max_discrete_velocity, min_mesh_size});
}

void BasicMesh::update_max_discrete_velocity(double value) {
    if (max_discrete_velocity <= 0.0) {
        max_discrete_velocity = value;
    } else {
        max_discrete_velocity = (max_discrete_velocity >= value) ? max_discrete_velocity : value;
    }
}

void BasicMesh::update_min_mesh_size(double value) {
    if (min_mesh_size <= 0.0) {
        min_mesh_size = value;
    } else {
        min_mesh_size = (min_mesh_size <= value) ? min_mesh_size : value;
    }
}


///     ListMesh
void ListMesh::load(const std::string &file_path) {
    su2Reader reader(file_path);
    if (!reader.is_file_open()) throw std::invalid_argument("open mesh file failed.");
    reader.parse(*this);
}

void ListMesh::build() {
    switch (type) {
        case MESH_TYPE_NORMAL:
            build_face();
        case MESH_TYPE_NO_FACE:
            build_geom();
        default:
            shrink_to_fit();
            {
                std::stringstream ss;
                ss << "Build mesh <" << name << "> - ok.";
                note_println(ss.str());
            }
            return;
    }
}

bool ListMesh::set_mark(const BoundaryParam &bc_param) {
    for (auto &mark : MARKS) {
        if (mark.name == bc_param.name) {
            // match
            mark.type = bc_param.type;
            mark.density = bc_param.density;
            mark.temperature = bc_param.temperature;
            mark.velocity = bc_param.velocity;
            for (auto &mark_elem : mark.MARK_ELEM) {
                auto &face = get_face(mark_elem.face_key);
                face.boundary_type = MESH::MarkTypeID.at(mark.type);
            }
            std::stringstream ss;
            ss << "Set mark params for Mesh<" << name << ">::Mark<" << mark.name << ">";
            info_println(ss.str());
            return true;
        }
    }
    // not match
    std::stringstream ss;
    ss << "Mark<" << bc_param.name << "> params have not been read from mesh file.";
    warn_println(ss.str());
    return false;
}

int ListMesh::node_num() {
    NPOIN = NODES.size();
    return NPOIN;
}

int ListMesh::cell_num() {
    NELEM = CELLS.size();
    return NELEM;
}

int ListMesh::face_num() const {
    return FACES.size();
}

Node<int> &ListMesh::get_node(const int &_key) {
    return NODES[_key];
}

Cell<int> &ListMesh::get_cell(const int &_key) {
    return CELLS[_key];
}

Face<int> &ListMesh::get_face(const int &_key) {
    return FACES[_key];
}

Mark<int> &ListMesh::get_mark(const int &_key) {
    return MARKS[_key];
}

void ListMesh::shrink_to_fit() {
    NODES.shrink_to_fit();
    CELLS.shrink_to_fit();
    FACES.shrink_to_fit();
    MARKS.shrink_to_fit();
    for (auto &cell : CELLS) cell.shrink_to_fit();
    for (auto &face : FACES) face.shrink_to_fit();
    for (auto &mark : MARKS) mark.shrink_to_fit();
}

void ListMesh::info() const {
    if (name == MESH_KEY_NULL) return;
    std::string type_name;
    switch (type) {
        case MESH_TYPE_NORMAL:
            type_name = "normal";
            break;
        case MESH_TYPE_NO_FACE:
            type_name = "no-face";
            break;
        default:
            type_name = "unknown";
    }
    std::stringstream ss;
    note_println("Mesh info:");
    ss << "  file: " << name << "\n"
       << "  type: " << type_name << "\n"
       << "  mesh data:";
    info_println(ss.str());
    data_int_println({"node", "cell", "face", "mark"}, {NPOIN, NELEM, int(FACES.size()), NMARK});
    data_double_println({"max-magnitude", "min-size"}, {max_discrete_velocity, min_mesh_size});
    if (type == MESH_TYPE_NORMAL) {
        info_println("  mark info:");
        for (auto &mark : MARKS) mark.info();
    }
}

///     MapMesh
void MapMesh::load(const std::string &file_path) {
    su2Reader reader(file_path);
    if (!reader.is_file_open()) throw std::invalid_argument("open mesh file failed.");
    reader.parse(*this);
}

void MapMesh::build() {
    switch (type) {
        case MESH_TYPE_NORMAL:
            build_face();
        case MESH_TYPE_NO_FACE:
            build_geom();
        default:
            shrink_to_fit();
            {
                std::stringstream ss;
                ss << "Build mesh <" << name << "> - ok.";
                note_println(ss.str());
            }
            return;
    }
}

bool MapMesh::set_mark(const BoundaryParam &bc_param) {
    for (auto &it : MARKS) {
        auto &mark_key = it.first;
        auto &mark = it.second;
        if (mark_key == bc_param.name) {
            // match
            mark.type = bc_param.type;
            mark.density = bc_param.density;
            mark.temperature = bc_param.temperature;
            mark.velocity = bc_param.velocity;
            for (auto &it2 : mark.MARK_ELEM) {
                auto &mark_elem = it2.second;
                auto &face = get_face(mark_elem.face_key);
                face.boundary_type = MESH::MarkTypeID.at(mark.type);
            }
            std::stringstream ss;
            ss << "Set mark params for Mark<" << mark.name << ">";
            info_println(ss.str());
            return true;
        }
    }
    // not match
    std::stringstream ss;
    ss << "Mark<" << bc_param.name << "> params have not been read from mesh file.";
    warn_println(ss.str());
    return false;
}

int MapMesh::node_num() {
    NPOIN = NodeKey.size();
    return NPOIN;
}

int MapMesh::cell_num() {
    NELEM = CellKey.size();
    return NELEM;
}

int MapMesh::face_num() const {
    return FaceKey.size();
}

Node<std::string> &MapMesh::get_node(const std::string &_key) {
    return NODES.at(_key);
}

Cell<std::string> &MapMesh::get_cell(const std::string &_key) {
    return CELLS.at(_key);
}

Face<std::string> &MapMesh::get_face(const std::string &_key) {
    return FACES.at(_key);
}

Mark<std::string> &MapMesh::get_mark(const std::string &_key) {
    return MARKS.at(_key);
}

void MapMesh::shrink_to_fit() {
    NodeKey.shrink_to_fit();
    CellKey.shrink_to_fit();
    FaceKey.shrink_to_fit();
    MarkKey.shrink_to_fit();
    /// rehash
    NODES.rehash(REHASH_SIZE);
    CELLS.rehash(REHASH_SIZE);
    FACES.rehash(REHASH_SIZE);
    MARKS.rehash(REHASH_SIZE);
    for (auto &it : CELLS) it.second.shrink_to_fit();
    for (auto &it : FACES) it.second.shrink_to_fit();
    for (auto &it : MARKS) it.second.shrink_to_fit();
}

void MapMesh::info() const {
    if (name == MESH_KEY_NULL) return;
    std::string type_name;
    switch (type) {
        case MESH_TYPE_NORMAL:
            type_name = "normal";
            break;
        case MESH_TYPE_NO_FACE:
            type_name = "no-face";
            break;
        default:
            type_name = "unknown";
    }
    std::stringstream ss;
    note_println("Mesh info:");
    ss << "  file: " << name << "\n"
       << "  type: " << type_name << "\n"
       << "  mesh data:";
    info_println(ss.str());
    data_int_println({"node", "cell", "face", "mark"}, {NPOIN, NELEM, int(FaceKey.size()), NMARK});
    data_double_println({"max-magnitude", "min-size"}, {max_discrete_velocity, min_mesh_size});
    if (type == MESH_TYPE_NORMAL) {
        info_println("  mark info:");
        for (auto &it : MARKS) it.second.info();
    }
}
