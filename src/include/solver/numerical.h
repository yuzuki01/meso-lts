/**
 * included by solver.h
 */

struct LeastSquare {
    std::vector<Vec3D> dr{};
    std::vector<double> weight{};
    Vec3D Cx{}, Cy{}, Cz{};
};

TP_key
LeastSquare generate_least_square(const key_type &cell_key, MESH::Mesh<key_type> &mesh);
