/**
 * included by solver.h
 */

/// LeastSquare
struct LeastSquare {
    std::vector<Vec3D> dr{};
    std::vector<double> weight{};
    // Vec3D Cx{}, Cy{}, Cz{};
    Mat3D C{};
};

TP_key_mesh
LeastSquare generate_least_square(const key_type &cell_key, mesh_type &mesh);

/// Ventaka limiter
struct Ventaka {
    double Wi, Wi_min, Wi_max, delta_ij;
};

double ventaka_phi(const Ventaka &_limiter, double omega);
