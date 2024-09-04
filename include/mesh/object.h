#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H

namespace MESO::fvmMesh {
    class Node;

    class Face;

    class Cell;

    class Mesh;

    struct LeastSquare;
    struct Symmetry;
}

namespace MESO {
    typedef std::vector<fvmMesh::Node> NodeList;
    typedef std::vector<fvmMesh::Face> FaceList;
    typedef std::vector<fvmMesh::Cell> CellList;

    struct fvmMesh::LeastSquare {
        /**
         * A * grad = b
         * grad = Inv(A) . b
         **/
        int neighbor_num{};
        ScalarList weight;
        VectorList dr;
        Vector Cx, Cy, Cz;
    };

    struct fvmMesh::Symmetry {
        int id[3]{};
    };
}

namespace MESO {

    enum {cell_field_flag, face_field_flag, node_field_flag};
    template<class FieldType>
    class Field {
    protected:
        fvmMesh::Mesh* mesh_ptr = nullptr;
    public:
        int len{};
        ObjectType flag{};
        typedef std::vector<FieldType> FieldArray;
        FieldArray values;

        /**
         * Realized in solver module.
         **/

        Field() = default;

        Field(fvmMesh::Mesh &mesh, int flag);

        Field(const Field<FieldType> &other);

        [[nodiscard]] fvmMesh::Mesh * get_mesh() const { return mesh_ptr; };

        Field<Scalar> heft(int id);

        FieldType &operator[](int index);

        void MeshCellValueToField(const std::function<FieldType(fvmMesh::Cell &)> &func);

        void set_zero();

        Field<Vector> gradient(bool _switch=true);

        void output(const std::string &file_name);
    };
}

namespace MESO::fvmMesh {
    class Node {
    public:
        ObjectId id;
        ObjectId group_id = 0;      /// 0 for interior node
        Position position;
        ObjectIdList neighbors;     /// neighbor cells contained the node

        Node(ObjectId id, const Position &position) : id(id), position(position) {};
    };

    class Face {
    public:
        ObjectType geom_type;
        ObjectId id;
        ObjectId group_id = 0;
        Position position;
        Scalar area{};
        ObjectIdList node_id;
        ObjectIdSet cell_id{-1, -1};    /// neighbor cell
        ObjectIdSet cell_face_id{-1, -1};
        VectorSet normal_vector;

        Face(ObjectId id, ObjectType geom_type, ObjectIdList node_list);
    };

    class Cell {
    public:
        ObjectType geom_type{};
        ObjectId id;
        ObjectId group_id = 0;
        ObjectId partition_id = 0;
        ObjectId partition_cell_id{};
        Position position;
        Scalar volume{};
        ObjectIdList node_id;
        ObjectIdList face_id;
        ObjectIdList neighbors;     /// neighbor cells shared the face
        LeastSquare least_square;
        Symmetry symmetry;

        explicit Cell(ObjectId id) : id(id) {};     // dvs init

        Cell(ObjectId id, ObjectType geom_type, ObjectIdList &node_list);

        void compute_least_square(const CellList &neighbor_cells, int dimension);  // defined in geom.cpp
    };

    class Mesh {
    public:
        int NNODE{}, NFACE{}, NCELL{}, NZONE{}, NMARK{}, NDFCD{}, NDFVL{};
        double total_volume{};
        NodeList nodes;
        FaceList faces;
        CellList cells;
        StringList cell_names;
        GroupList cell_groups;
        StringList face_names = {"fluid-interior"};
        GroupList face_groups;
        GroupList cell_partition_groups, face_partition_groups;

        double min_cell_size{};
        double max_cell_magnitude{};

        Mesh() = default;

        [[nodiscard]] int dimension() const { return NDFVL; };

        void generate_face();

        void update_mesh_params();

        void build_geom(double scale_ratio=1.0);

#ifdef _METIS_H_
        void partition();
#endif

        void output(const std::string &file_name);

        void output(const std::string &file_name,
                    std::initializer_list<std::string> names,
                    std::initializer_list<Field<Scalar> *> values,
                    int step=-1,double solution_time=-1.0);

        Field<Scalar> zero_scalar_field(int flag=cell_field_flag);
        Field<Vector> zero_vector_field(int flag=cell_field_flag);

        void info();
    };
}

template<class FieldType>
void MESO::Field<FieldType>::MeshCellValueToField(const std::function<FieldType(fvmMesh::Cell &)> &func) {
    for (int i = 0; i < len; ++i) {
        values[i] = func(mesh_ptr->cells[i]);
    }
}

#endif //MESH_OBJECT_H
