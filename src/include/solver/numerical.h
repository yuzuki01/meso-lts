/**
 * included by solver.h
 */

/// LeastSquare
struct LeastSquare {
    std::vector<Vec3D> dr{};
    std::vector<double> weight{};
    Mat3D C;
};

TP_mesh
LeastSquare generate_least_square(int cell_key, mesh_type &mesh);

/// Ventaka limiter
struct Ventaka {
    double Wi, Wi_min, Wi_max, delta_ij;
};

double ventaka_phi(const Ventaka &_limiter, double omega);
