/***********************************
 *            Math Lib             *
 ***********************************
 *  1. 3D Vector <-                *
 *  2. 3 x 3 Matrix                *
 ***********************************/

#include "core.h"

/// 3D Vec3D

double Vec3D::operator[](int i) const {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("Index out of range for Vec3D.");
    }
}

double &Vec3D::operator[](int i) {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("Index out of range for Vec3D.");
    }
}

Vec3D &Vec3D::operator=(const Vec3D &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

Vec3D &Vec3D::operator=(const std::initializer_list<double> &init_list) {
    switch (init_list.size()) {
        case 3:
            z = *(init_list.begin() + 2);
        case 2:
            y = *(init_list.begin() + 1);
        case 1:
            x = *init_list.begin();
            break;
        default:
            throw std::invalid_argument("Vec3D: initializer_list has a invalid size.");
    }
    return *this;
}

Vec3D Vec3D::operator+(const Vec3D &other) const {
    // 相加
    return Vec3D(x + other.x, y + other.y, z + other.z);
}

Vec3D Vec3D::operator-(const Vec3D &other) const {
    // 相减
    return Vec3D(x - other.x, y - other.y, z - other.z);
}

Vec3D Vec3D::operator-() const {
    // 取反
    return Vec3D(-x, -y, -z);
}

Vec3D operator*(double k, const Vec3D &v) {
    // 乘于系数
    return Vec3D(v.x * k, v.y * k, v.z * k);
}

double operator*(const Vec3D &v1, const Vec3D &v2) {
    // 点乘
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec3D Vec3D::operator*(double k) const {
    // 乘于系数
    return Vec3D(x * k, y * k, z * k);
}

Vec3D Vec3D::operator/(double k) const {
    // 除于系数
    if (k != 0.0) {
        return Vec3D(x / k, y / k, z / k);
    } else {
        return Vec3D(INFINITY, INFINITY, INFINITY);
    }
}

Vec3D &Vec3D::operator+=(Vec3D &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec3D &Vec3D::operator-=(Vec3D &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec3D &Vec3D::operator+=(const Vec3D &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec3D &Vec3D::operator-=(const Vec3D &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec3D &Vec3D::operator*=(double k) {
    x *= k;
    y *= k;
    z *= k;
    return *this;
}

Vec3D &Vec3D::operator/=(double k) {
    x /= k;
    y /= k;
    z /= k;
    return *this;
}

double Vec3D::operator*(Vec3D &other) const {
    // 点乘
    return x * other.x + y * other.y + z * other.z;
}

Vec3D Vec3D::operator^(Vec3D &other) const {
    // 实现叉乘运算
    double result_x = y * other.z - z * other.y;
    double result_y = z * other.x - x * other.z;
    double result_z = x * other.y - y * other.x;
    return Vec3D(result_x, result_y, result_z);
}

double Vec3D::magnitude() const {
    // 取模
    return sqrt(x * x + y * y + z * z);
}

double Vec3D::sum() const {
    return x + y + z;
}

Vec3D Vec3D::norm() const {
    // 归一
    double k = magnitude();
    // 除于系数
    return *this / k;
}

std::string Vec3D::info() const {
    std::stringstream ss;
    ss << "<" << std::setprecision(OUT_PRECISION) << x << ", " << y << ", " << z << ">";
    return ss.str();
}
