/**
 * included by mesh.h
 */

#define MESH_KEY_NULL "NULL"

#define MESH_TYPE_NORMAL    0
#define MESH_TYPE_NO_FACE   1

/// template
#define TP_func template<>
#define TP_key template<class key_type>
#define TP_mesh template<class mesh_type>
#define TP_key_mesh template<class key_type, class mesh_type>
#define key_vector std::vector<key_type>

namespace MESH {

    TP_key
    class Node;

    TP_key
    class Cell;

    TP_key
    class Face;

    TP_key
    class MarkElem;

    TP_key
    class Mark;

    class BasicMesh;    // basic mesh class
    /// Mesh
    class ListMesh;
    class MapMesh;

    extern std::unordered_map<std::string, int> MarkTypeID;
}

struct BoundaryParam {
    std::string name{};
    std::string type{};
    double density = 0.0;
    double temperature = 0.0;
    Vec3D velocity = {0.0, 0.0, 0.0};
};

/// MESH geom object
TP_key
class MESH::Node {
public:
    key_type key{};
    Vec3D position = {0.0, 0.0, 0.0};

    explicit Node(const key_type &node_key, const Vec3D &node_position);
};

TP_key
class MESH::Cell {
public:
    int type;           /// Geom type
    key_type key;

    /// for structural DVS - bounce-back
    key_type inv_cell_key{};
    // Cell<key_type> * inv_cell_ptr;
    /// for structural DVS - mirror-reflect
    int on_layer_id{}, layer_cell_num{};
    /// Geom parameters
    int face_num;                       /// number of interface
    bool have_shadow;                   /// shadow cell
    key_type shadow_cell_key{};         /// key of shadow cell
    Vec3D position = {0.0, 0.0, 0.0};
    double position_square;             /// equal to position * position
    double volume;                      /// volume of cell
    key_vector node_key;
    // std::vector<Node<key_type> *> node_ptr;
    key_vector face_key;
    // std::vector<Face<key_type> *> face_ptr;
    key_vector near_cell_key;           /// cells using the same face
    // std::vector<Cell<key_type> *> near_cell_ptr;
    key_vector second_near_cell_key;    /// for high-order scheme
    // std::vector<Cell<key_type> *> second_near_cell_ptr;

    void shrink_to_fit();

    /// su2 format              TYPE | NODES | (ID)
    explicit Cell(const key_type &cell_key, const int cell_type, const key_vector &node_set);

    /// structural DVS mesh
    Cell(const Vec3D &particle_velocity, double weight,
         const key_type &cell_key, const key_type &inv_cell_key,
         int on_layer_cell_id, int on_layer_cell_num);

    /// debug func
    std::string info() const;

    std::string info_with_pos() const;
};

TP_key
class MESH::Face {
public:
    int type{};           /// Geom type
    key_type key{};
    int boundary_type = MESH_BC_INTERFACE;  /// boundary type
    key_type mark_key{};                    /// mark where the interface is
    // Mark<key_type> * mark_ptr{};

    /// Geom parameters
    /**
     * on_cell_key -> cell_1
     *  on_cell_nv -> normal Vec3D pointing from face center to cell-1 center.
     *  on_cell_face -> cell_1.face_key[on_cell_face] == this->key
     * inv_cell_key -> cell_2
     *  inv_cell_nv -> normal Vec3D pointing from face center to cell_2 center.
     *  inv_cell_face -> cell_2.face_key[inv_cell_face] == this->key
     *  When the interface is on boundary (boundary type == 0):
     *      on_cell_key == inv_cell_key
     *      on_cell_nv == inv_cell_nv
     */
    key_type on_cell_key{};
    // Cell<key_type> * on_cell_ptr{};
    int on_cell_face{};
    Vec3D on_cell_nv{};
    key_type inv_cell_key{};
    // Cell<key_type> * inv_cell_ptr{};
    int inv_cell_face{};
    Vec3D inv_cell_nv{};
    double area{};
    Vec3D position = {0.0, 0.0, 0.0};
    key_vector node_key{};
    // std::vector<Node<key_type> *> node_ptr{};

    void shrink_to_fit();

    /// self defined format       face_id | face_type | node_id
    Face(const key_type &face_key, int elem_type, const key_vector &node_key_vec);

    Face() = default;

    /// debug func
    std::string info() const;

    std::string info_all() const;
};

TP_key
class MESH::MarkElem {
public:
    int type;
    std::string key;
    key_vector node_set;
    key_type face_key{};

    /// su2 format              TYPE | nodes
    explicit MarkElem(const string_vector &init_strings);
};

TP_func
class MESH::Mark<int> {
public:
    std::string name;
    int elem_num;
    std::string type = MESH_KEY_NULL;
    /// Physical parameters
    double density = 0.0, temperature = 0.0;
    Vec3D velocity = {0.0, 0.0, 0.0};

    /// su2
    explicit Mark(const std::string &mark_name, const int nelem);

    void shrink_to_fit();

    /// get item
    MarkElem<int> &get_elem(const int &_key);

    /// debug func
    void info() const;

    std::vector<MarkElem<int>> MARK_ELEM;
};

TP_func
class MESH::Mark<std::string> {
public:
    std::string name;
    int elem_num;
    std::string type = MESH_KEY_NULL;
    /// Physical parameters
    double density = 0.0, temperature = 0.0;
    Vec3D velocity = {0.0, 0.0, 0.0};

    /// su2
    explicit Mark(const std::string &mark_name, const int nelem);

    void shrink_to_fit();

    /// get item
    MarkElem<std::string> &get_elem(const std::string &_key);

    /// debug func
    void info() const;

    std::vector<std::string> MarkElemKey;
    std::unordered_map<std::string, MarkElem<std::string>> MARK_ELEM;
};

/// MESH basic class
class MESH::BasicMesh {
public:
    /**
     * MapMesh params
     */
    int type = -1;
    int NDIME = 0;  // dimension
    int NPOIN = 0;  // number of node points
    int NELEM = 0;  // number of elements
    int NMARK = 0;  // number of marks
    std::string name{};
    double max_discrete_velocity = -1.0;
    double min_mesh_size = -1.0;

    BasicMesh() = default;

    explicit BasicMesh(int mesh_type) : type(mesh_type) {};

    BasicMesh(int mesh_type, std::string mesh_name) : type(mesh_type), name(std::move(mesh_name)) {};

    int dimension() const;

    void update_max_discrete_velocity(double _value);

    void update_min_mesh_size(double _value);

    void info();
};

/// ListMesh
class MESH::ListMesh : public BasicMesh {
public:
    /**
     * Container
     */
    std::vector<Node<int>> NODES{};
    std::vector<Cell<int>> CELLS{};
    std::vector<Face<int>> FACES{};
    std::vector<Mark<int>> MARKS{};

    ListMesh() = default;

    explicit ListMesh(int mesh_type) : BasicMesh(mesh_type) {};

    ListMesh(int mesh_type, std::string mesh_name) : BasicMesh(mesh_type, mesh_name) {};

    void load(const std::string &file_path);

    void build();

    void set_mark(const BoundaryParam &bc_param);

    void info() const;

    void shrink_to_fit();

    /// get item
    int node_num();

    int cell_num();

    int face_num() const;

    Node<int> &get_node(const int &_key);

    Cell<int> &get_cell(const int &_key);

    Face<int> &get_face(const int &_key);

    Mark<int> &get_mark(const int &_key);

private:
    /// mesh/mesh_build.cpp
    void build_face();       // create interface
    void build_geom();            // calculate geom params
};

/// MapMesh
class MESH::MapMesh : public BasicMesh {
public:
    /**
     * Container
     */
    /// key
    std::vector<std::string> NodeKey{};
    std::vector<std::string> CellKey{};
    std::vector<std::string> FaceKey{};
    std::vector<std::string> MarkKey{};
    /// map
    std::unordered_map<std::string, Node<std::string>> NODES{};
    std::unordered_map<std::string, Cell<std::string>> CELLS{};
    std::unordered_map<std::string, Face<std::string>> FACES{};
    std::unordered_map<std::string, Mark<std::string>> MARKS{};

    MapMesh() = default;

    explicit MapMesh(int mesh_type) : BasicMesh(mesh_type) {};

    MapMesh(int mesh_type, std::string mesh_name) : BasicMesh(mesh_type, mesh_name) {};

    void load(const std::string &file_path);

    void build();

    void set_mark(const BoundaryParam &bc_param);

    void info() const;

    void shrink_to_fit();

    /// get item
    int node_num();

    int cell_num();

    int face_num() const;

    Node<std::string> &get_node(const std::string &_key);

    Cell<std::string> &get_cell(const std::string &_key);

    Face<std::string> &get_face(const std::string &_key);

    Mark<std::string> &get_mark(const std::string &_key);

private:
    /// mesh/mesh_build.cpp
    void build_face();       // create interface
    void build_geom();            // calculate geom params
};
