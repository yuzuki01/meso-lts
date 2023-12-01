#include "mesh.h"

/// local function
TP_mesh
bool is_in_cell_tria(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    /// basic shape
    /// lp = u * l1 + v * l2
    double u, v, C;
    auto &n0 = mesh.get_node(node_set[0]).position,
            &n1 = mesh.get_node(node_set[1]).position,
            &n2 = mesh.get_node(node_set[2]).position;
    Vec3D l1 = n1 - n0, l2 = n2 - n0, lp = p_center - n0;
    C = l1.y * l2.x - l1.x * l2.y;
    u = (l2.x * lp.y - l2.y * lp.x) / C;
    v = (l1.y * lp.x - l1.x * lp.y) / C;
    return u > 0.0 && v > 0.0 && (u + v) < 1.0;
}

TP_mesh
bool is_in_cell_quad(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    return is_in_cell_tria<mesh_type>(p_center, {node_set[0], node_set[1], node_set[2]}, mesh) ||
           is_in_cell_tria<mesh_type>(p_center, {node_set[0], node_set[3], node_set[2]}, mesh);
}

TP_mesh
bool is_in_cell_hexah(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    return false;
}

TP_mesh
bool is_in_cell_pyram(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    return false;
}

TP_mesh
bool is_in_cell_tetra(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    return false;
}

TP_mesh
bool is_in_cell_prism(const Vec3D &p_center, const std::vector<int> &node_set, mesh_type &mesh) {
    return false;
}

TP_mesh
bool is_particle_in_cell(const Vec3D &p_center, MESH::Cell &cell, mesh_type &mesh) {
    std::vector<int> node_set(GEOM::node_num(cell.type));
    switch (cell.type) {
        case GEOM::TRIA:
            return is_in_cell_tria(p_center, node_set, mesh);
        case GEOM::QUAD:
            return is_in_cell_quad(p_center, node_set, mesh);
        case GEOM::TETRA:
            return is_in_cell_tetra(p_center, node_set, mesh);
        case GEOM::HEXAH:
            return is_in_cell_hexah(p_center, node_set, mesh);
        case GEOM::PRISM:
            return is_in_cell_prism(p_center, node_set, mesh);
        case GEOM::PYRAM:
            return is_in_cell_pyram(p_center, node_set, mesh);
        default:
            return false;
    }
}

/// DPM::Particle
TP_func bool DPM::Particle<MESH::StaticMesh>::is_in_cell(MESH::Cell &cell) {
    return is_particle_in_cell(position, cell, mesh);
}

TP_func bool DPM::Particle<MESH::MapMesh>::is_in_cell(MESH::Cell &cell) {
    return is_particle_in_cell(position, cell, mesh);
}
