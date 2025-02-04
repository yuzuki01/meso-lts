#ifndef CORE_MATH_H
#define CORE_MATH_H

#define VSMALL 1e-30

namespace MESO::Math {
    class Vector;
}

class MESO::Math::Vector {
public:
    double x, y, z;
    Vector() : x(0.0), y(0.0), z(0.0) {};
    Vector(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
    Vector(const std::initializer_list<double> &init_list);
    Vector& operator=(const Vector &other);
    double& operator[](int index);
    double operator[](int index) const;

    Vector& operator+=(const Vector& other);
    Vector& operator-=(const Vector& other);
    Vector& operator*=(double _s);
    Vector& operator/=(double _s);

    Vector operator-() const;
    Vector operator*(double _s) const;
    Vector operator/(double _s) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    double operator*(const Vector& other) const;
    Vector operator^(const Vector& other) const;

    [[nodiscard]] double magnitude() const;
    [[nodiscard]] Vector normalize() const;
    [[nodiscard]] std::string str() const;
};

MESO::Math::Vector operator*(double k, const MESO::Math::Vector &vec);

/// openMP
#pragma omp declare reduction(+: MESO::Math::Vector : omp_out += omp_in)
#pragma omp declare reduction(-: MESO::Math::Vector : omp_out -= omp_in)

#endif //CORE_MATH_H
