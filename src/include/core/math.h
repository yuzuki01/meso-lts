#ifndef CORE_MATH_H
#define CORE_MATH_H

#define VSMALL 1e-30

namespace MESO::Math {
    class Vector;
}

class MESO::Math::Vector {
public:
    Scalar x, y, z;
    Vector() : x(0.0), y(0.0), z(0.0) {};
    explicit Vector(Scalar _v) : x(_v), y(_v), z(_v) {};
    Vector(Scalar _x, Scalar _y, Scalar _z) : x(_x), y(_y), z(_z) {};
    Vector(const std::initializer_list<Scalar> &init_list);
    Vector& operator=(const Vector &other);
    Scalar& operator[](int index);
    Scalar operator[](int index) const;

    Vector& operator+=(const Vector& other);
    Vector& operator-=(const Vector& other);
    Vector& operator*=(Scalar _s);
    Vector& operator/=(Scalar _s);

    Vector operator-() const;
    Vector operator*(Scalar _s) const;
    Vector operator/(Scalar _s) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Scalar operator*(const Vector& other) const;
    Vector operator^(const Vector& other) const;

    [[nodiscard]] Scalar magnitude() const;
    [[nodiscard]] Vector normalize() const;
    [[nodiscard]] std::string str() const;
};

MESO::Math::Vector operator*(MESO::Scalar k, const MESO::Math::Vector &vec);

MESO::Scalar mag(const MESO::Math::Vector &x);

MESO::Scalar magSqr(const MESO::Math::Vector &x);

/// openMP
#pragma omp declare reduction(+: MESO::Math::Vector : omp_out += omp_in)
#pragma omp declare reduction(-: MESO::Math::Vector : omp_out -= omp_in)

#endif //CORE_MATH_H
