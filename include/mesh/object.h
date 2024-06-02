#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H

namespace MESO::Mesh {
    class Node;

    class Face;

    class Cell;

    class Zone;

    struct LeastSquare;
}

namespace MESO {
    typedef std::vector<Mesh::Node> NodeList;
    typedef std::vector<Mesh::Face> FaceList;
    typedef std::vector<Mesh::Cell> CellList;

    struct Mesh::LeastSquare {
        /**
         * A * grad = b
         * grad = Inv(A) . b
         **/
        int neighbor_num{};
        ScalarList weight;
        VectorList dr;
        Vector Cx, Cy, Cz;
    };
}

namespace MESO {

    enum {cell_field_flag, face_field_flag};
    template<class FieldType>
    class Field {
    protected:
        Mesh::Zone* mesh_ptr = nullptr;
    public:
        int len{};
        ObjectType flag{};
        typedef std::vector<FieldType> FieldArray;
        FieldArray values;

        /**
         * Realized in solver module.
         **/

        Field() = default;

        Field(Mesh::Zone &mesh, int flag);

        Field(const Field<FieldType> &other);

        [[nodiscard]] Mesh::Zone * get_mesh() const { return mesh_ptr; };

        Field<Scalar> heft(int id);

        FieldType &operator[](int index);

        void MeshCellValueToField(const std::function<FieldType(Mesh::Cell &)> &func);

        Field<Vector> gradient(bool _switch=true);
    };
}

namespace MESO::Mesh {
    class Node {
    public:
        ObjectId id;
        Position position;

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
        ObjectIdSet on_cell_id{-1, -1};
        VectorSet normal_vector;

        Face(ObjectId id, ObjectType geom_type, ObjectIdList node_list);
    };

    class Cell {
    public:
        ObjectType geom_type{};
        ObjectId id;
        ObjectId group_id = 0;
        Position position;
        Scalar volume{};
        ObjectIdList node_id;
        ObjectIdList face_id;
        ObjectIdList neighbors;     // neighbor cells shared faces
        LeastSquare least_square;

        explicit Cell(ObjectId id) : id(id) {};     // dvs init

        Cell(ObjectId id, ObjectType geom_type, ObjectIdList &node_list);

        void calculate_least_square(const CellList &neighbor_cells, int dimension);  // defined in geom.cpp
    };

    class Zone {
    public:
        int NNODE{}, NFACE{}, NCELL{}, NZONE{}, NMARK{}, NDFCD{}, NDFVL{};
        NodeList nodes;
        FaceList faces;
        CellList cells;
        StringList cell_names;
        GroupList cell_groups;
        StringList face_names = {"fluid-interior"};
        GroupList face_groups;

        double min_cell_size{};
        double max_cell_magnitude{};

        Zone() = default;

        [[nodiscard]] int dimension() const { return NDFVL; };

        void generate_face();

        void update_num();

        void build_geom();

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

#endif //MESH_OBJECT_H
