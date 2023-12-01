/***********************************
 *            Math Lib             *
 ***********************************
 *  1. 3D Vector                   *
 *  2. 3 x 3 Matrix  <-            *
 ***********************************/

#include "core.h"

Mat3D::Mat3D(double every) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            _matrix[i][j] = every;
        }
    }
}

Mat3D::Mat3D(std::initializer_list<double> values) {
    if (values.size() != 9)
        warn_println("Warning: The provided initializer list does not have length 9. "
                     "The matrix will be partially filled, and the excess values will be ignored or set to 0.");

    auto it = values.begin();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (it != values.end()) {
                _matrix[i][j] = *it;
                ++it;
            } else {
                _matrix[i][j] = 0.0;
            }
        }
    }
}

double Mat3D::det() const {
    return _matrix[0][0] * (_matrix[1][1] * _matrix[2][2] - _matrix[2][1] * _matrix[1][2]) -
           _matrix[0][1] * (_matrix[1][0] * _matrix[2][2] - _matrix[2][0] * _matrix[1][2]) +
           _matrix[0][2] * (_matrix[1][0] * _matrix[2][1] - _matrix[2][0] * _matrix[1][1]);
}

Mat3D Mat3D::T() {
    Mat3D result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result._matrix[i][j] = _matrix[j][i];
        }
    }
    return result;
}

Mat3D Mat3D::I() {
    Mat3D result;
    // 计算矩阵的行列式
    double _det = det();
    if (_det == 0.0) {
        warn_println("Error: The matrix is singular, and its inverse does not exist.");
        // 处理矩阵不可逆的情况
        // 这里你可以添加适当的处理逻辑，比如返回一个特殊值表示不可逆，或者抛出异常等。
    } else {
        // 计算伴随矩阵
        result._matrix[0][0] = (_matrix[1][1] * _matrix[2][2] - _matrix[2][1] * _matrix[1][2]) / _det;
        result._matrix[0][1] = -(_matrix[0][1] * _matrix[2][2] - _matrix[2][1] * _matrix[0][2]) / _det;
        result._matrix[0][2] = (_matrix[0][1] * _matrix[1][2] - _matrix[1][1] * _matrix[0][2]) / _det;
        result._matrix[1][0] = -(_matrix[1][0] * _matrix[2][2] - _matrix[2][0] * _matrix[1][2]) / _det;
        result._matrix[1][1] = (_matrix[0][0] * _matrix[2][2] - _matrix[2][0] * _matrix[0][2]) / _det;
        result._matrix[1][2] = -(_matrix[0][0] * _matrix[1][2] - _matrix[1][0] * _matrix[0][2]) / _det;
        result._matrix[2][0] = (_matrix[1][0] * _matrix[2][1] - _matrix[2][0] * _matrix[1][1]) / _det;
        result._matrix[2][1] = -(_matrix[0][0] * _matrix[2][1] - _matrix[2][0] * _matrix[0][1]) / _det;
        result._matrix[2][2] = (_matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1]) / _det;
    }
    return result;
}

Mat3D Mat3D::operator*(const Mat3D &_mat) {
    Mat3D result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result._matrix[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                result._matrix[i][j] += _matrix[i][k] * _mat._matrix[k][j];
            }
        }
    }
    return result;
}

Vec3D Mat3D::operator*(const Vec3D &_vec) {
    double x = 0.0, y = 0.0, z = 0.0;
    for (int j = 0; j < 3; ++j) {
        x += _matrix[0][j] * _vec[j];
        y += _matrix[1][j] * _vec[j];
        z += _matrix[2][j] * _vec[j];
    }
    return Vec3D{x, y, z};
}

void Mat3D::info() const {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << std::setw(15) << std::scientific << _matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
}
