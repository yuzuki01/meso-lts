#include "core.h"

Mat3D rotate_matrix(const Vec3D &vec) {
    /**
     * Input:   Vec3D   vec
     * Output:  Mat3D   M
     *      where M satisfies the equations:
     *          n = vec.norm()
     *          M * n = {1,0,0}
     *          M.I() * {1,0,0} = n
     **/
    Vec3D n = vec.norm();
    Vec3D v = {1.0, 0.0, 0.0};
    Vec3D k = n ^ v;
    Mat3D M;
    double k_len = k.magnitude();
    if (k_len == 0.0) {
        M = identity();
        if (n * v < 0.0) M = -M;
    } else {
        k = k / k_len;
        double theta = acos(n * v);
        Mat3D I = identity();
        Mat3D kk{k.x * k.x, k.x * k.y, k.x * k.z,
                 k.x * k.y, k.y * k.y, k.y * k.z,
                 k.x * k.z, k.y * k.z, k.z * k.z,};
        Mat3D K{0, -k.z, k.y,
                k.z, 0, -k.x,
                -k.y, k.x, 0};
        M = cos(theta) * I + (1.0 - cos(theta)) * kk + sin(theta) * K;
    }
    return M;
}
