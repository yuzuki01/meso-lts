# meso/core

meso 核心

## 参数解析器 Argparser
```c++
int main() {
    /// 实例化解析器
    Argparser parse(argc, argv);
    
    /// 解析启动参数 --max_step <value> 解析失败则返回 "100"
    auto it1 = parser.parse_param("max_step", "100");
    int max_step = stoi(it1);
    /// 解析启动参数 -debug 有则返回 true，无则返回 false
    auto it2 = parser.parse_switch("debug");
    bool debug_mode = it2;
    
    /// other codes
    
    return 0;
}
```

## 日志对象 Logger
```c++
int main() {
    /// 日志初始化
    LoggerInit();
    /// 实例化日志对象
    Logger logger("LoggerName");
    
    /**
     * 流式输出
     * 结果：
     *  $ [LoggerName] This is a test line.
     */
    logger << "This is a test line.";
    logger.debug();
    
    /// other codes
    
    return 0;
}
```

## 强制输出
```c++
info_println("Hello, world!");
note_println("Hello, world!");
```

## 自定义数学库 mesomath.h

### Vec3D

```c++
/// 三维向量
///     vec1 = <1, 2, 3>
///     vec2 = <0, 0, 0> (default)
Vec3D vec1(1.0, 2.0, 3.0), vec2;
/// 向量相加减
Vec3D plus = vec1 + vec2;
Vec3D ngtv = vec1 - vec2;
/// 向量放缩
double scale = 2.0;
Vec3D multi = scale * vec1;
Vec3D devid = vec2 / scale;
/// 向量点乘
double dot = vec1 * vec2;
/// 向量叉乘
Vec3D cross = vec1 ^ vec2;
/// 向量自加自减
vec1 += vec2;
vec2 -= vec1;
/// 正则化
vec1.norm();
/// 模长
double magnitude = vec2.magnitude();
```

### Mat3D

```c++
/// 3x3矩阵
/// mat1 = {
//  1, 0, 0,
//  0, 1, 0,
//  0, 0, 1,
// };
Mat3D mat1{1,0,0,0,1,0,0,0,1};
double det = mat1.det();    /// 计算行列式
Mat3D tra_mat = mat1.T();    /// 转置
Mat3D inv_mat = mat1.I();    /// 求逆
Vec3D vec1{1,2,3};, vec2;
vec2 = mat * vec1;          /// 矩阵乘于列向量
Mat3D mat2{3,0,0,0,2,0,0,0,1};
Mat3D mulpti_mat = mat1 * mat2; /// 矩阵 mat1 乘于 mat2
```

## 工具 utils

### split

```c++
/// .h
#define string_vector std::vector<std::string>

/// .cpp
std::string line = "a b  c";
string_vector result = split(line);
/**
 * Output:
 * $ a
 * $ b
 * $ c
 */
for (auto &it : result) {
    std::cout << it << std::endl;
}
```

### str_vec_cmp
```c++
string_vector line = {"a", "b", "c"};

str_vec_cmp(line, "a b c d");   /// false
str_vec_cmp(line, "b b c");     /// false
str_vec_cmp(line, "a d c");     /// false
str_vec_cmp(line, "a b c");     /// true
str_vec_cmp(line, "a b");       /// true
str_vec_cmp(line, "a");         /// true
```
