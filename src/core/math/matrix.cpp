#include "core/core.h"

using namespace MESO::Math;


Matrix::Matrix() {
    for (auto &vec: data) {
        vec = {0.0, 0.0, 0.0};
    }
}

Matrix::Matrix(const std::initializer_list<double> &init_list) {
    const int len = int(init_list.size());
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int index = i * 3 + j;
            data[i][j] = (index < len) ? *(init_list.begin() + index) : 0.0;
        }
    }
}

Vector &Matrix::operator[](int row) {
    return data[row];
}

Matrix &Matrix::operator+=(MESO::Math::Matrix &other) {
    for (int i = 0; i < 3; ++i) {
        data[i] += other[i];
    }
    return *this;
}


Matrix &Matrix::operator-=(MESO::Math::Matrix &other) {
    for (int i = 0; i < 3; ++i) {
        data[i] -= other[i];
    }
    return *this;
}

Matrix &Matrix::operator*=(double _s) {
    for (auto &vec: data) {
        vec *= _s;
    }
    return *this;
}

Matrix &Matrix::operator/=(double _s) {
    double _ss = 1.0 / _s;
    for (auto &vec: data) {
        vec *= _ss;
    }
    return *this;
}

Vector Matrix::operator*(const Vector& vec) {
    return {data[0] * vec, data[1] * vec, data[2] * vec};
}

Matrix Matrix::operator*(double _s) {
    Matrix result(*this);
    for (auto &vec: result.data) {
        vec *= _s;
    }
    return result;
}

Matrix Matrix::operator/(double _s) {
    Matrix result(*this);
    double _ss = 1.0 / _s;
    for (auto &vec: result.data) {
        vec *= _ss;
    }
    return result;
}

Matrix Matrix::operator+(MESO::Math::Matrix &other) {
    Matrix result(*this);
    for (int i = 0; i < 3; ++i) {
        result[i] += other[i];
    }
    return result;
}

Matrix Matrix::operator-(MESO::Math::Matrix &other) {
    Matrix result(*this);
    for (int i = 0; i < 3; ++i) {
        result[i] -= other[i];
    }
    return result;
}

Matrix Matrix::operator*(const MESO::Math::Matrix &other) {
    Matrix result, other_t = other.T();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = data[i] * other_t[j];
        }
    }
    return result;
}

double Matrix::Det() const {
    auto &v1 = data[0];
    auto &v2 = data[1];
    auto &v3 = data[2];
    return -v1.z * v2.y * v3.x + v1.y * v2.z * v3.x
           + v1.z * v2.x * v3.y - v1.x * v2.z * v3.y
           - v1.y * v2.x * v3.z + v1.x * v2.y * v3.z;
}

Matrix Matrix::T() const {
    Matrix result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = data[j][i];
        }
    }
    return result;
}

Matrix Matrix::I() const {
    Matrix result;
    auto &v1 = data[0];
    auto &v2 = data[1];
    auto &v3 = data[2];
    double fm = -v1.z * v2.y * v3.x + v1.y * v2.z * v3.x
                + v1.z * v2.x * v3.y - v1.x * v2.z * v3.y
                - v1.y * v2.x * v3.z + v1.x * v2.y * v3.z;
    result[0] = {v2.y * v3.z - v2.z * v3.y, v1.z * v3.y - v1.y * v3.z, v1.y * v2.z - v1.z * v2.y};
    result[1] = {v2.z * v3.x - v2.x * v3.z, v1.x * v3.z - v1.z * v3.x, v1.z * v2.x - v1.x * v2.z};
    result[2] = {v2.x * v3.y - v2.y * v3.x, v1.y * v3.x - v1.x * v3.y, v1.x * v2.y - v1.y * v2.x};
    return result / fm;
}

MESO::String Matrix::str() {
    std::stringstream ss;
    ss << "\n|" << std::setw(7) << std::left << data[0][0] << " " << std::setw(7) << std::left << data[0][1] << " "
       << std::setw(7) << std::left << data[0][2] << "|\n"
       << "|" << std::setw(7) << std::left << data[1][0] << " " << std::setw(7) << std::left << data[1][1] << " "
       << std::setw(7) << std::left << data[1][2] << "|\n"
       << "|" << std::setw(7) << std::left << data[2][0] << " " << std::setw(7) << std::left << data[2][1] << " "
       << std::setw(7) << std::left << data[2][2] << "|";
    return ss.str();
}

MESO::Math::Matrix operator*(double k, MESO::Math::Matrix mat) {
    return mat * k;
}
