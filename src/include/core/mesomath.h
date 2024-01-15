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

    double operator[](int i) const;
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
    double sum() const;                             // 求和
    Vec3D norm() const;                                   // 归一

    std::string info() const;
};

Vec3D operator*(double k, const Vec3D &v);            // 乘于系数
double operator*(const Vec3D &v1, const Vec3D &v2);   // 点乘

/// 3D Matrix
class Mat3D {
public:
    double _matrix[3][3]{0.0};
    ~Mat3D() = default;
    Mat3D() = default;
    explicit Mat3D(double every);    // 矩阵各个值设置为 every
    explicit Mat3D(std::initializer_list<double> values);
    Mat3D &operator=(const Mat3D &_mat);

    double det() const;
    Mat3D T();                              // 转置
    Mat3D I();                              // 矩阵的逆
    Mat3D operator*(double _value);         // 乘于系数
    Mat3D operator/(double _value);         // 除于系数
    Vec3D operator*(const Vec3D &_vec);     // Mat * Vec 矩阵乘于列向量
    Mat3D operator*(const Mat3D &_mat);     // 矩阵乘法 this * _mat
    Mat3D operator+(const Mat3D &_mat);
    Mat3D operator-();
    Mat3D operator-(const Mat3D &_mat);

    void info() const;
};

Mat3D operator*(double _value, const Mat3D &_mat);
Mat3D identity();
Mat3D diag(std::initializer_list<double> values);

/// 3D Coordinate Transfer
Mat3D rotate_matrix(const Vec3D &_vec);
