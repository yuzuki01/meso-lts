#include "mesh.h"

/// local function
TP_key_mesh
bool is_in_cell_tria(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
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

TP_key_mesh
bool is_in_cell_quad(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
    return is_in_cell_tria<key_type, mesh_type>(p_center, {node_set[0], node_set[1], node_set[2]}, mesh) ||
           is_in_cell_tria<key_type, mesh_type>(p_center, {node_set[0], node_set[3], node_set[2]}, mesh);
}

TP_key_mesh
bool is_in_cell_hexah(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
    return false;
}

TP_key_mesh
bool is_in_cell_pyram(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
    return false;
}

TP_key_mesh
bool is_in_cell_tetra(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
    return false;
}

TP_key_mesh
bool is_in_cell_prism(const Vec3D &p_center, const std::vector<key_type> &node_set, mesh_type &mesh) {
    return false;
}

TP_key_mesh
bool is_particle_in_cell(const Vec3D &p_center, MESH::Cell<key_type> &cell, mesh_type &mesh) {
    std::vector<key_type> node_set(GEOM::node_num(cell.type));
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
TP_func bool DPM::Particle<int, MESH::ListMesh>::is_in_cell(MESH::Cell<int> &cell) {
    return is_particle_in_cell(position, cell, mesh);
}

TP_func bool DPM::Particle<std::string, MESH::MapMesh>::is_in_cell(MESH::Cell<std::string> &cell) {
    return is_particle_in_cell(position, cell, mesh);
}
