#ifndef CORE_MATH_H
#define CORE_MATH_H

namespace MESO::Math {
    class Vector;
    class Matrix;
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

class MESO::Math::Matrix {
public:
    std::array<Vector, 3> data;
    Matrix();
    Matrix(const Matrix &other) = default;
    Matrix& operator=(const Matrix& other) = default;
    Matrix(const std::initializer_list<double> &init_list);
    Vector& operator[](int row);
    Matrix& operator+=(Matrix& other);
    Matrix& operator-=(Matrix& other);
    Matrix& operator*=(double _s);
    Matrix& operator/=(double _s);

    Vector operator*(const Vector& vec);
    Matrix operator*(double _s);
    Matrix operator/(double _s);

    Matrix operator+(Matrix& other);
    Matrix operator-(Matrix& other);
    Matrix operator*(const Matrix& other);

    [[nodiscard]] double Det() const;
    [[nodiscard]] Matrix T() const;
    [[nodiscard]] Matrix I() const;
    std::string str();
};

MESO::Math::Matrix operator*(double k, MESO::Math::Matrix mat);

#endif //CORE_MATH_H
