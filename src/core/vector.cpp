#include "core/core.h"
#include "core/math.h"


namespaceMESO
using namespace Math;


Vector::Vector(const std::initializer_list<double> &init_list) {
    const int len = int(init_list.size());
    if (len >= 3) z = *(init_list.begin() + 2); else z = 0.0;
    if (len >= 2) y = *(init_list.begin() + 1); else y = 0.0;
    if (len >= 1) x = *init_list.begin(); else x = 0.0;
}

Vector &Vector::operator=(const Math::Vector &other) = default;


double &Vector::operator[](int index) {
    switch (index) {
        case X_COMPONENT:
            return x;
        case Y_COMPONENT:
            return y;
        case Z_COMPONENT:
            return z;
        default:
            throw std::invalid_argument("<Math::Vector> index out of range.");
    }
}

double Vector::operator[](int index) const {
    switch (index) {
        case X_COMPONENT:
            return x;
        case Y_COMPONENT:
            return y;
        case Z_COMPONENT:
            return z;
        default:
            throw std::invalid_argument("<Math::Vector> index out of range.");
    }
}

Vector &Vector::operator+=(const Math::Vector &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector &Vector::operator-=(const Math::Vector &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vector &Vector::operator*=(double _s) {
    x *= _s;
    y *= _s;
    z *= _s;
    return *this;
}

Vector &Vector::operator/=(double _s) {
    double _ss = 1.0 / _s;
    x *= _ss;
    y *= _ss;
    z *= _ss;
    return *this;
}

Vector Vector::operator-() const {
    return {-x, -y, -z};
}

Vector Vector::operator*(double _s) const {
    return {x * _s, y * _s, z * _s};
}

Vector Vector::operator/(double _s) const {
    double _ss = 1.0 / _s;
    return {x * _ss, y * _ss, z * _ss};
}

Vector Vector::operator+(const Math::Vector &other) const {
    return {x + other.x, y + other.y, z + other.z};
}

Vector Vector::operator-(const Math::Vector &other) const {
    return {x - other.x, y - other.y, z - other.z};
}

double Vector::operator*(const Vector &other) const {
    return x * other.x + y * other.y + z * other.z;
}

Vector Vector::operator^(const Math::Vector &other) const {
    return {y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x};
}

double Vector::magnitude() const {
    return sqrt(x * x + y * y + z * z);
}

Vector Vector::normalize() const {
    double mag = this->magnitude();
    if (mag == 0.0) {
        logger.warn << "<Math::Vector> Normalized a zero-vector." << std::endl;
        return {0.0, 0.0, 0.0};
    }
    return *this / mag;
}

std::string Vector::str() const {
    std::stringstream ss;
    ss << "[" << x << "," << y << "," << z << "]";
    return ss.str();
}

Math::Vector operator*(Scalar k, const Math::Vector &vec) {
    return {vec.x * k, vec.y * k, vec.z * k};
}


Scalar mag(const Math::Vector &x) {
    return x.magnitude();
}


Scalar magSqr(const Math::Vector &x) {
    return x * x;
}
