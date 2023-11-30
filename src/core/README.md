# meso/core

meso 核心

## 参数解析器 Argparser

参数解析对象，解析两种形式的启动参数：

 - 开关类型的参数，格式为 `-<switch_name>` ，启动参数中含有则返回 `true` 否则为 `false` ；
 - 键值类型的参数，格式为 `--<key> <value>` ，对应函数需要传入 value 类型和默认值。
   
示例代码如下：

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

流形式的日志输出对象，示例代码如下：

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

带有日志等级颜色格式的输出，该函数会在结尾默认附带一个 `\n`，传入的字符串不需要额外添加换行符。

```c++
info_println("Hello, world!");
note_println("Hello, world!");
```

## 自定义数学库 mesomath.h

### Vec3D

三维向量类，支持数学上的矢量加减，矢量放缩，矢量点乘和矢量叉乘操作。

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

数学上的 3x3 矩阵，支持行列式运算，转置，求逆。

与 `Vec3D` 相乘时，按照 `mat * vec` 的顺序，`vec` 被视为一个数学上的列向量，结果得到一个列向量。

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

将 ```std::string``` 按照分隔符 `\n`  `\t` 或空格拆分为 ```std::vector<std::string>``` ，主要用于 `Reader` 类的参数解析。

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

该函数传入参数为 `str_vec_cmp(const string_vector &_vec, const std::string &reference)` ，其作用是将传入的 `_vec` 与通过 `split(std::string)` 函数处理的 `std::string reference` 进行比较。

比较规则为，首先以 `std::string reference` 分割后得到的 `string_vector` 类型的结果为基准。当传入的 `_vec` 的长度小于基准的长度时，返回 `false` ；当长度相当或大于时，采取逐项对比的方法，超出基准长度的部分被忽略，若被比较的部分逐项均相等则返回 `true` 。代码示例如下：

```c++
string_vector line = {"a", "b", "c"};

str_vec_cmp(line, "a b c d");   /// false
str_vec_cmp(line, "b b c");     /// false
str_vec_cmp(line, "a d c");     /// false
str_vec_cmp(line, "a b c");     /// true
str_vec_cmp(line, "a b");       /// true
str_vec_cmp(line, "a");         /// true
```
