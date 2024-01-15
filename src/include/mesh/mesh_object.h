#include <utility>

/**
 * included by mesh.h
 */

#define MESH_KEY_NULL "NULL"

// #define MESH_TYPE_NORMAL    0
// #define MESH_TYPE_NO_FACE   1
enum {MeshTypeNormal, MeshTypeNoFace};

/// template
#define TP_func template<>
#define TP_mesh template<class mesh_type>
#define key_vector std::vector<int>

namespace MESH {

    class Node;

    class Cell;

    class Face;

    class MarkElem;

    class Mark;

    class BasicMesh;    // basic mesh class
    /// Mesh
    class StaticMesh;
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
class MESH::Node {
public:
    int key{};
    Vec3D position = {0.0, 0.0, 0.0};

    explicit Node(int node_key, const Vec3D &node_position);
};

class MESH::Cell {
public:
    int type;           /// Geom type
    int key;

    /// for structural DVS - bounce-back
    int inv_cell_key{};
    // Cell<int> * inv_cell_ptr;
    /// for structural DVS - mirror-reflect
    int on_layer_id{}, layer_cell_num{};
    /// Geom parameters
    int face_num;                           /// number of interfaces
    int ghost_cell_num = 0;                 /// number of ghost cells
    int near_ghost_cell_key{};
    key_vector second_near_ghost_cell_key{};        /// key of second-near ghost cells
    Vec3D position = {0.0, 0.0, 0.0};
    double position_square;             /// equal to position * position
    double volume;                      /// volume of cell
    key_vector node_key;
    // std::vector<Node<int> *> node_ptr;
    key_vector face_key;
    // std::vector<Face<int> *> face_ptr;
    key_vector near_cell_key;           /// cells using the same face
    // std::vector<Cell<int> *> near_cell_ptr;
    key_vector second_near_cell_key;    /// for high-order scheme
    // std::vector<Cell<int> *> second_near_cell_ptr;

    void shrink_to_fit();

    /// su2 format              TYPE | NODES | (ID)
    explicit Cell(int cell_key, int cell_type, const key_vector &node_set);

    /// structural DVS mesh
    Cell(const Vec3D &particle_velocity, double weight,
         const int &cell_key, const int &inv_cell_key,
         int on_layer_cell_id, int on_layer_cell_num);

    /// debug func
    std::string info() const;

    std::string info_with_pos() const;
};

class MESH::Face {
public:
    int type{};           /// Geom type
    int key{};
    int boundary_type = MeshBC_interface;  /// boundary type
    int mark_key{};                    /// mark where the interface is
    // Mark<int> * mark_ptr{};

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
    int on_cell_key{};
    // Cell<int> * on_cell_ptr{};
    int on_cell_face{};
    Vec3D on_cell_nv{};
    Mat3D rotate_matrix{}, inv_rotate_matrix{};     /// for on_cell only
    int inv_cell_key{};
    // Cell<int> * inv_cell_ptr{};
    int inv_cell_face{};
    Vec3D inv_cell_nv{};
    double area{};
    Vec3D position = {0.0, 0.0, 0.0};
    key_vector node_key{};
    // std::vector<Node<int> *> node_ptr{};

    void shrink_to_fit();

    /// self defined format       face_id | face_type | node_id
    Face(int face_key, int elem_type, const key_vector &node_key_vec);

    Face() = default;

    /// debug func
    std::string info() const;

    std::string info_all() const;
};

class MESH::MarkElem {
public:
    int type;
    std::string key;
    key_vector node_set;
    int face_key{};

    /// su2 format              TYPE | nodes
    explicit MarkElem(const string_vector &init_strings);
};

class MESH::Mark {
public:
    std::string name;
    int elem_num;
    std::string type = MESH_KEY_NULL;
    /// Physical parameters
    double density = 0.0, temperature = 0.0;
    Vec3D velocity = {0.0, 0.0, 0.0};

    /// su2
    explicit Mark(const std::string &mark_name, int nelem);

    void shrink_to_fit();

    /// get item
    MarkElem &get_elem(int _key);

    /// debug func
    void info() const;

    std::vector<MarkElem> MARK_ELEM;
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

/// StaticMesh
class MESH::StaticMesh : public BasicMesh {
public:
    /**
     * Container
     */
    std::vector<Node> NODES{};
    std::vector<Cell> CELLS{};
    std::vector<Face> FACES{};
    std::vector<Mark> MARKS{};

    StaticMesh() = default;

    explicit StaticMesh(int mesh_type) : BasicMesh(mesh_type) {};

    StaticMesh(int mesh_type, std::string mesh_name) : BasicMesh(mesh_type, std::move(mesh_name)) {};

    void load(const std::string &file_path);

    void build();

    bool set_mark(const BoundaryParam &bc_param);

    void compute_rotate_matrix();

    void info() const;

    void shrink_to_fit();

    /// get item
    int node_num();

    int cell_num();

    int face_num() const;

    Node &get_node(int _key);

    Cell &get_cell(int _key);

    Face &get_face(int _key);

    Mark &get_mark(int _key);

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
    std::vector<int> NodeKey{};
    std::vector<int> CellKey{};
    std::vector<int> FaceKey{};
    std::vector<int> MarkKey{};
    /// map
    std::unordered_map<int, Node> NODES{};
    std::unordered_map<int, Cell> CELLS{};
    std::unordered_map<int, Face> FACES{};
    std::unordered_map<int, Mark> MARKS{};

    MapMesh() = default;

    explicit MapMesh(int mesh_type) : BasicMesh(mesh_type) {};

    MapMesh(int mesh_type, std::string mesh_name) : BasicMesh(mesh_type, std::move(mesh_name)) {};

    void load(const std::string &file_path);

    void build();

    bool set_mark(const BoundaryParam &bc_param);

    void info() const;

    void shrink_to_fit();

    /// get item
    int node_num();

    int cell_num();

    int face_num() const;

    Node &get_node(int _key);

    Cell &get_cell(int _key);

    Face &get_face(int _key);

    Mark &get_mark(int _key);

private:
    /// mesh/mesh_build.cpp
    void build_face();       // create interface
    void build_geom();       // calculate geom params
};
