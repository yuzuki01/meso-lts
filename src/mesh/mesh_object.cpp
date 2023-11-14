#include "mesh.h"

using namespace MESH;

/// global vars
std::unordered_map<std::string, int> MESH::MarkTypeID {
        {"interface", MESH_BC_INTERFACE},
        {"inlet", MESH_BC_INLET},
        {"outlet", MESH_BC_OUTLET},
        {"isothermal_wall", MESH_BC_ISOTHERMAL_WALL},
        {"adiabat_wall", MESH_BC_ABDIABAT_WALL},
        {"symmetry", MESH_BC_SYMMETRY},
        {"periodic", MESH_BC_PERIODIC},
};

/// Node
TP_func Node<int>::Node(const std::vector<std::string> &init_strings) {
    const int len = init_strings.size();
    if (len >= 3) {
        key = stoi(init_strings[0]) - 1;
        position.x = stod(init_strings[1]);
        position.y = stod(init_strings[2]);
        if (len >= 4) {
            position.z = stod(init_strings[3]);
        } else {
            position.z = 0.0;
        }
    } else {
        throw std::invalid_argument("neu_node_formatter: invalid init_strings.size().");
    }
}

TP_func Node<std::string>::Node(const std::vector<std::string> &init_strings) {
    const int len = init_strings.size();
    if (len >= 3) {
        key = init_strings[0];
        position.x = stod(init_strings[1]);
        position.y = stod(init_strings[2]);
        if (len >= 4) {
            position.z = stod(init_strings[3]);
        } else {
            position.z = 0.0;
        }
    } else {
        throw std::invalid_argument("neu_node_formatter: invalid init_strings.size().");
    }
}

/// Cell
TP_func Cell<int>::Cell(const std::vector<std::string> &init_strings) {
    key = stoi(init_strings[0]) - 1;
    type = stoi(init_strings[1]);
    const int node_num = stoi(init_strings[2]);
    for (int i = 0; i < node_num; i++) {
        node_key.push_back(stoi(init_strings[i + 3]) - 1);
    }
    if (int(node_key.size()) != GEOM::node_num(type)) {
        std::stringstream ss;
        ss << "Cell<key=" << key << "> caught invalid node key set.";
        error_println(ss.str());
        throw std::invalid_argument("Cell caught invalid node key set.");
    }
    volume = -1.0;
    have_shadow = false;
    face_num = GEOM::face_num(type);
    face_key.resize(face_num, -1);
    position_square = -1.0;
}

TP_func Cell<std::string>::Cell(const std::vector<std::string> &init_strings) {
    key = init_strings[0];
    type = stoi(init_strings[1]);
    const int node_num = stoi(init_strings[2]);
    for (int i = 0; i < node_num; i++) {
        node_key.push_back(init_strings[i + 3]);
    }
    if (int(node_key.size()) != GEOM::node_num(type)) {
        std::stringstream ss;
        ss << "Cell<key=" << key << "> caught invalid node key set.";
        error_println(ss.str());
        throw std::invalid_argument("Cell caught invalid node key set.");
    }
    volume = -1.0;
    have_shadow = false;
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
    have_shadow = false;
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
    have_shadow = false;
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
    cell_key = stoi(init_strings[0]) - 1;
    cell_type = stoi(init_strings[1]);
    /// NOTE:   .neu index starts from 1
    ///         std::vector index starts from 0
    on_cell_face_id = stoi(init_strings[2]) - 1;
}

TP_func MarkElem<std::string>::MarkElem(const std::vector<std::string> &init_strings) {
    cell_key = init_strings[0];
    cell_type = stoi(init_strings[1]);
    /// NOTE:   .neu index starts from 1
    ///         std::vector index starts from 0
    on_cell_face_id = stoi(init_strings[2]) - 1;
}

/// Mark
Mark<int>::Mark(const string_vector &init_strings) {
    const int len = init_strings.size();
    if (len > 0) key = init_strings[0];
    if (len > 2) elem_num = stoi(init_strings[2]);
}

Mark<std::string>::Mark(const string_vector &init_strings) {
    const int len = init_strings.size();
    if (len > 0) key = init_strings[0];
    if (len > 2) elem_num = stoi(init_strings[2]);
}

void Mark<int>::shrink_to_fit() {
    MARK_ELEM.shrink_to_fit();
}

void Mark<std::string>::shrink_to_fit() {
    MarkElemKey.shrink_to_fit();
    MARK_ELEM.rehash(REHASH_SIZE);
}

MarkElem<int> & Mark<int>::get_elem(const int &_key) {
    return MARK_ELEM[_key];
}

MarkElem<std::string> & Mark<std::string>::get_elem(const std::string &_key) {
    return MARK_ELEM.at(_key);
}

void Mark<int>::info() const {
    std::stringstream ss;
    ss << "  Mark<" << key <<">\n"
       << "     - elements: " << elem_num << "\n"
       << "     - type: " << type << "\n"
       << "     - density: " << density << " temperature: " << temperature << "\n"
       << "     - velocity: " << velocity.info();
    info_println(ss.str());
}

void Mark<std::string>::info() const {
    std::stringstream ss;
    ss << "  Mark<" << key <<">\n"
       << "     - elements: " << elem_num << "\n"
       << "     - type: " << type << "\n"
       << "     - density: " << density << " temperature: " << temperature << "\n"
       << "     - velocity: " << velocity.info();
    info_println(ss.str());
}

/// Mesh
///     Basic mesh:
int BasicMesh::dimension() const {
    return NDFCD;
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
    data_int_println({"node", "cell", "face", "mark"}, {NUMNP, NELEM, -1, NBSETS});
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
    neuReader reader(file_path);
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

void ListMesh::set_mesh_params(const std::vector<std::string> &init_strings) {
    NUMNP = stoi(init_strings[0]);
    NELEM = stoi(init_strings[1]);
    NGRPS = stoi(init_strings[2]);
    NBSETS = stoi(init_strings[3]);
    NDFCD = stoi(init_strings[4]);
    NDFVL = stoi(init_strings[5]);
    // reserve
    NODES.reserve(NUMNP);
    CELLS.reserve(NELEM);
    MARKS.reserve(NBSETS);
}

void ListMesh::set_mark(const BoundaryParam &bc_param) {
    for (auto &mark : MARKS) {
        if (mark.key == bc_param.name) {
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
            ss << "Set mark params for Mark<" << mark.key << ">";
            info_println(ss.str());
            return;
        }
    }
    // not match
    std::stringstream ss;
    ss << "Mark<" << bc_param.name << "> params have not been read from mesh file.";
    warn_println(ss.str());
}

int ListMesh::node_num() {
    NUMNP = NODES.size();
    return NUMNP;
}

int ListMesh::cell_num() {
    NELEM = CELLS.size();
    return NELEM;
}

int ListMesh::face_num() const {
    return FACES.size();
}

Node<int> & ListMesh::get_node(const int &_key) {
    return NODES[_key];
}

Cell<int> & ListMesh::get_cell(const int &_key) {
    return CELLS[_key];
}

Face<int> & ListMesh::get_face(const int &_key) {
    return FACES[_key];
}

Mark<int> & ListMesh::get_mark(const int &_key) {
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
    data_int_println({"node", "cell", "face", "mark"}, {NUMNP, NELEM, int(FACES.size()), NBSETS});
    data_double_println({"max-magnitude", "min-size"}, {max_discrete_velocity, min_mesh_size});
    if (type == MESH_TYPE_NORMAL) {
        info_println("  mark info:");
        for (auto &mark : MARKS) mark.info();
    }
}

///     MapMesh
void MapMesh::load(const std::string &file_path) {
    neuReader reader(file_path);
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

void MapMesh::set_mesh_params(const std::vector<std::string> &init_strings) {
    NUMNP = stoi(init_strings[0]);
    NELEM = stoi(init_strings[1]);
    NGRPS = stoi(init_strings[2]);
    NBSETS = stoi(init_strings[3]);
    NDFCD = stoi(init_strings[4]);
    NDFVL = stoi(init_strings[5]);
    // reserve
    NodeKey.reserve(NUMNP);
    CellKey.reserve(NELEM);
    MarkKey.reserve(NBSETS);
    NODES.reserve(NUMNP);
    CELLS.reserve(NELEM);
    MARKS.reserve(NBSETS);
}

void MapMesh::set_mark(const BoundaryParam &bc_param) {
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
            ss << "Set mark params for Mark<" << mark.key << ">";
            info_println(ss.str());
            return;
        }
    }
    // not match
    std::stringstream ss;
    ss << "Mark<" << bc_param.name << "> params have not been read from mesh file.";
    warn_println(ss.str());
}

int MapMesh::node_num() {
    NUMNP = NodeKey.size();
    return NUMNP;
}

int MapMesh::cell_num() {
    NELEM = CellKey.size();
    return NELEM;
}

int MapMesh::face_num() const {
    return FaceKey.size();
}

Node<std::string> & MapMesh::get_node(const std::string &_key) {
    return NODES.at(_key);
}

Cell<std::string> & MapMesh::get_cell(const std::string &_key) {
    return CELLS.at(_key);
}

Face<std::string> & MapMesh::get_face(const std::string &_key) {
    return FACES.at(_key);
}

Mark<std::string> & MapMesh::get_mark(const std::string &_key) {
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
    data_int_println({"node", "cell", "face", "mark"}, {NUMNP, NELEM, int(FaceKey.size()), NBSETS});
    data_double_println({"max-magnitude", "min-size"}, {max_discrete_velocity, min_mesh_size});
    if (type == MESH_TYPE_NORMAL) {
        info_println("  mark info:");
        for (auto &it : MARKS) it.second.info();
    }
}
