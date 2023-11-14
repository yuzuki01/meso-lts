/**
 * included by core.h
 *  Math Lib
 */

/// 3D Vector
class Vec3D {
public:
    double x, y, z;

    Vec3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    Vec3D() : x(0.0), y(0.0), z(0.0) {}

    ~Vec3D() = default;

    double &operator[](int i);                      // 赋值
    Vec3D &operator=(const Vec3D &other);

    Vec3D &operator=(const std::initializer_list<double> &init_list);

    Vec3D operator+(const Vec3D &other) const;      // 相加
    Vec3D operator-(const Vec3D &other) const;      // 相减
    Vec3D operator-() const;                        // 取反
    double operator*(Vec3D &other) const;           // 点乘
    Vec3D operator^(Vec3D &other) const;            // 叉乘
    Vec3D operator*(double k) const;                // 乘于系数
    Vec3D operator/(double k) const;                // 除于系数
    Vec3D &operator+=(Vec3D &other);                // 自加
    Vec3D &operator-=(Vec3D &other);                // 自减
    Vec3D &operator+=(const Vec3D &other);          // 自加
    Vec3D &operator-=(const Vec3D &other);          // 自减
    Vec3D &operator*=(double k);                    // 自乘 double
    Vec3D &operator/=(double k);                    // 自除
    double magnitude() const;                       // 取模
    void norm();                                    // 归一

    std::string info() const;
};

Vec3D operator*(double k, const Vec3D &v);            // 乘于系数
double operator*(const Vec3D &v1, const Vec3D &v2);   // 点乘
