#ifndef MESH_GEOM_H
#define MESH_GEOM_H

namespace MESO::Geom {
    const ObjectType Edge = 1;
    const ObjectType Quad = 2;
    const ObjectType Tria = 3;
    const ObjectType Brick = 4;
    const ObjectType Wedge = 5;
    const ObjectType Tetra = 6;
    const ObjectType Pyram = 7;

    int face_num(ObjectType geom_type);

    int node_num(ObjectType geom_type);

    Vector calculate_position(const List<fvmMesh::Node> &node_list);

    Scalar calculate_area(ObjectType geom_type, const List<fvmMesh::Node> &node_list);

    Scalar calculate_volume(ObjectType geom_type, const List<fvmMesh::Node> &node_list);

    Vector calculate_normal_vector(ObjectType face_geom_type, Vector &cell_position, Vector &face_position,
                                   const List<fvmMesh::Node> &node_list);

    KeyString generate_key(const List<ObjectId> &node_list);
}

#endif //MESH_GEOM_H
