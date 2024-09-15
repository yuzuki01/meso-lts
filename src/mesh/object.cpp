#include "mesh/mesh.h"


using namespace MESO::fvmMesh;


Face::Face(MESO::ObjectId id, MESO::ObjectType geom_type, MESO::List<ObjectId> node_list)
        : id(id), geom_type(geom_type), node_id(std::move(node_list)) {}

Cell::Cell(MESO::ObjectId id, MESO::ObjectType geom_type, MESO::List<ObjectId> &node_list)
        : id(id), geom_type(geom_type), node_id(std::move(node_list)) {
    face_id.resize(Geom::face_num(geom_type), -1);
}

/**
 * generate face
 **/
using namespace MESO::Geom;

MESO::Map<MESO::List<MESO::ObjectId>> fg_map = {
        {Quad,  {Edge, Edge, Edge, Edge}},
        {Tria,  {Edge, Edge, Edge}},
        {Brick, {Quad, Quad, Quad, Quad, Quad, Quad}},
        {Wedge, {Quad, Quad, Quad, Tria, Tria}},
        {Tetra, {Tria, Tria, Tria, Tria}},
        {Pyram, {Quad, Tria, Tria, Tria, Tria}}
};

MESO::Map<MESO::List<MESO::List<MESO::ObjectId>>> fni_map = {
        {Quad,  {{0, 1},       {1, 2},       {2, 3},       {3, 0}}},
        {Tria,  {{0, 1},       {1, 2},       {2, 0}}},
        {Brick, {{0, 1, 5, 4}, {1, 3, 7, 5}, {3, 2, 6, 7}, {2, 0, 4, 6}, {1, 0, 2, 3}, {4, 5, 7, 6,}}},
        {Wedge, {{0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {0, 2, 1},    {3, 4, 5}}},
        {Tetra, {{1, 0, 2},    {0, 1, 3},    {1, 2, 3},    {2, 0, 3}}},
        {Pyram, {{0, 2, 3, 1}, {0, 1, 4},    {1, 3, 4},    {3, 2, 4},    {2, 0, 4}}}
};

void Mesh::generate_face() {
    Dict<ObjectId> fm_edge, fm_quad, fm_tria;
    for (auto &cell: cells) {
        for (int fi = 0; fi < face_num(cell.geom_type); ++fi) {
            /// range faces on cell
            ObjectType geom_type = fg_map[cell.geom_type][fi];
            List<ObjectId> node_list;
            for (auto it: fni_map[cell.geom_type][fi]) {
                node_list.push_back(cell.node_id[it]);
            }
            String face_key = generate_key(node_list);
            Dict<ObjectId> *fmp;
            switch (geom_type) {
                case Edge:
                    fmp = &fm_edge;
                    break;
                case Quad:
                    fmp = &fm_quad;
                    break;
                case Tria:
                    fmp = &fm_tria;
                    break;
                default:
                    throw std::invalid_argument("<generate_face> unsupported face.");
            }
            Dict<ObjectId> &fm = *fmp;
            /// check face
            auto it = fm.find(face_key);
            if (it == fm.end()) {
                /// not exist
                int face_id = int(faces.size());
                faces.emplace_back(face_id, geom_type, node_list);
                fm[face_key] = face_id;
                auto &face = faces.back();
                face.cell_id = {cell.id, cell.id};
                face.cell_face_id = {fi, fi};
                cell.face_id[fi] = face.id;
            } else {
                auto &face = faces[it->second];
                face.cell_id[1] = cell.id;
                face.cell_face_id[1] = fi;
                cell.face_id[fi] = face.id;
            }
        }
    }
    NFACE = int(faces.size());
}

void Mesh::update_mesh_params() {
    NNODE = int(nodes.size());
    NFACE = int(faces.size());
    NCELL = int(cells.size());
    NZONE = int(cell_groups.size());
    NMARK = int(face_groups.size());
}

void Mesh::build_geom(double scale_ratio) {
    /// scale
    for (auto &node : nodes) {
        node.position *= scale_ratio;
    }
    /// update mesh params
    update_mesh_params();
    /// cell
    total_volume = 0.0;
    const int D = dimension();
    min_cell_size = -1.0;
    max_cell_magnitude = -1.0;
    cell_partition_groups.resize(MPI::processor_num, {});
    for (auto &cell: cells) {
        cell.partition_id = MPI::rank;
        cell.partition_cell_id = cell.id;
        cell_partition_groups[MPI::rank].push_back(cell.id);
        List<Node> node_list;
        for (auto &it: cell.node_id) {
            node_list.push_back(nodes[it]);
        }
        cell.position = Geom::calculate_position(node_list);
        cell.volume = Geom::calculate_volume(cell.geom_type, node_list);
        total_volume += cell.volume;
        double mcs = (D == 2) ? 2.0 * sqrt(cell.volume / M_PI) : 2.0 *
                                                                 pow((3.0 * cell.volume / (4.0 * M_PI)), 1.0 / 3.0);
        double mcm = cell.position.magnitude();
        if (min_cell_size < 0.0) {
            min_cell_size = mcs;
        } else {
            min_cell_size = (mcs < min_cell_size) ? mcs : min_cell_size;
        }
        if (max_cell_magnitude < 0.0) {
            max_cell_magnitude = mcm;
        } else {
            max_cell_magnitude = (mcm > max_cell_magnitude) ? mcm : max_cell_magnitude;
        }
    }
    cell_partition_groups.shrink_to_fit();
    /// face
    for (auto &face: faces) {
        List<Node> node_list;
        for (auto &it: face.node_id) {
            auto &node = nodes[it];
            node_list.push_back(node);
            node.group_id = face.group_id;
        }
        face.position = Geom::calculate_position(node_list);
        face.area = Geom::calculate_area(face.geom_type, node_list);
        face.normal_vector[0] = Geom::calculate_normal_vector(face.geom_type, cells[face.cell_id[0]].position,
                                                              face.position, node_list);
        face.normal_vector[1] = -face.normal_vector[0];
    }
    /// neighbor
    for (auto &cell: cells) {
        List<Cell> neighbors;
        for (int face_id: cell.face_id) {
            auto &face = faces[face_id];
            for (int cell_id: face.cell_id) {
                if (cell_id == cell.id) continue;
                cell.neighbors.push_back(cell_id);
                neighbors.push_back(cells[cell_id]);
            }
        }
        /// least square
        cell.compute_least_square(neighbors, dimension());
    }
    /// shrink_to_fit
    nodes.shrink_to_fit();
    faces.shrink_to_fit();
    cells.shrink_to_fit();
    cell_names.shrink_to_fit();
    cell_groups.shrink_to_fit();
    face_names.shrink_to_fit();
    face_groups.shrink_to_fit();
}

void Mesh::info() {
    logger.note << "fvmMesh info: \n";
    logger.info << "    Node num: " << NNODE << "\n"
                << "    Face num: " << NFACE << "\n"
                << "    Cell num: " << NCELL << "\n"
                << "    Min-cell-size: " << min_cell_size << "\n"
                << "    Max-cell-magnitude: " << max_cell_magnitude << "\n";
    logger.info << "    Cell Groups: \n";
    for (int i = 0; i < NZONE; ++i) {
        logger.info << "        " << i + 1 << " " << cell_names[i] << "\n"
                    << "            cell num: " << int(cell_groups[i].size()) << "\n";
    }
    logger.info << "    Face Groups: \n";
    for (int i = 0; i < NMARK; ++i) {
        logger.info << "        " << i + 1 << " " << face_names[i] << "\n"
                    << "            face num: " << int(face_groups[i].size()) << "\n";
    }
    logger.info << std::endl;
}
