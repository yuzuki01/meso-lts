# meso/mesh

meso 网格

---

## 目录
 
 - [网格分类](#网格分类)
 - [宏定义](#宏定义)
 - [网格对象](#网格对象)
    - [网格结点](#网格结点)
    - [网格单元](#网格单元)
    - [网格界面](#网格界面)
 - [网格生成器](#网格生成器)
    - [Gauss-Hermit型DVS](#Gauss-Hermit型DVS)
    - [Newton-Cotes型DVS](#Newton-Cotes型DVS)

---

## 网格分类

当前网格分类：

- 固定网格 `MESH::StaticMesh` ，其容器为 `std::vector<MESH::OBJECT>`

- 可变网格 `MESH::MapMesh` ，其容器为 `std::unordered_map<int, MESH::OBJECT>`

不同的网格类通过模板函数实现了通用接口

---

## 宏定义

定义了三个与网格模板相关的宏，其中 `int` 表示网格对象的 `key` 值的类型。

`int` 的类型为固定网格 `int` ，可变网格 `std::string` 。

```c++
#define TP_mesh template <class mesh_type>
#define std::vector<int> key_vector
#define std::vector<std::string> string_vector
```

---

## 网格对象

### 网格结点

```c++
/// 网格结点
class MESH::Node {
public:
    // key 值
    int key;
    // 坐标
    Vec3D position = {0.0, 0.0, 0.0};
    // 按照 cpp 结构构造
    explicit Node(int node_key, const Vec3D &node_position);
    // 按照 su2 格式构造
    explicit Node(const string_vector &init_strings);
};
```

### 网格单元

```c++
/// 网格单元
class MESH::Cell {
public:
    // su2 几何体类型
    int type;
    // key 值
    int key;
    // 对于结构速度空间 - 反方向 cell 的 key 值
    int inv_cell_key;
    // 对于结构速度空间 - 镜面反射
    int on_layer, layer_num;
    // 几何参数
    // 界面个数
    int face_num;
    // 影单元
    bool have_shadow;
    int shadow_cell_key;
    // 坐标
    Vec3D position;
    // 对于速度空间 坐标模长平方
    double position_square;
    // 体积
    double volume;
    // 构成单元的节点 key 向量
    key_vector node_key;
    // 构成单元的界面 key 向量
    key_vector face_key;
    // 相邻单元的 key 向量
    key_vector near_cell_key;
    // 次相邻单元的 key 相量    for high-order scheme
    key_vector second_near_cell_key;
    /// 按照 su2 格式构造
    explicit Cell(const string_vector &init_strings);
    /// 结构速度空间构造函数
    Cell(const Vec3D &particle_velocity, double weight,
         int cell_key, const int &inv_cell_key,
         int on_layer_cell_id, int on_layer_cell_num);
};
```

### 网格界面

```c++
/// 网格界面
class MESH::Face {
public:
    // su2 几何体类型
    int type;
    // key 值
    int key;
    // 边界类型     默认为 0 ，非边界界面
    int boundary_type = 0;
    // 界面所在边界(MESH::Mark) key 值
    std::string mark_key = MESH_KEY_NULL;

    /// 几何参数
    /**
     * on_cell 表示正向 inv_cell 表示反向，正反向的判断仅取决于界面构建时的顺序，无特殊含义。
     *  界面为边界时，由于只有一个相邻单元(MESH::Cell)，
     *  on_cell 对应的值绝对代表了界面所附着的唯一的网格单元
     *  且此时有如下的一些变量关系：
     *      on_cell_key == inv_cell_key         反向单元就是正向单元
     *      on_cell_face == inv_cell_face       反向单元就是正向单元
     *      on_cell_nv == -inv_cell_nv          法向量依旧存在取反关系
     *      on_cell_nv 为单元格心向外的界面法向量
     *      inv_cell_nv 为边界指向流体域的界面法向量
     */
    
    int on_cell_key;   // 正向单元的 key 值
    int on_cell_face;       // 正向单元上的界面编号
    Vec3D on_cell_nv;       // 正向单元指向外的界面法向量
    int inv_cell_key;  // 反向单元的 key 值
    int inv_cell_face;      // 反向单元上的界面编号
    Vec3D inv_cell_nv;      // 反向单元指向外的界面法向量
    // 面积
    double area;
    // 坐标
    Vec3D position = {0.0, 0.0, 0.0};
    // 构成界面的节点(MESH::Node) key 值
    key_vector node_key;

    /// 构造
    Face(int face_key, int elem_type, const key_vector &node_key_vec);
    Face();
    
    // get item
    // 模板函数统一调用接口
    Node &get_node(int _key);

    Cell &get_cell(int _key);

    Face &get_face(int _key);

    Mark &get_mark(int _key);
};
```

---

## 网格生成器

### Gauss-Hermit型DVS

生成 Gauss-Hermit 型结构速度空间的函数。其中 `gauss_point` 按照一维情况来给，函数会自动根据网格文件的维数来生成对应的结构型速度空间网格。

如： `dvs_mesh = GENERATOR::gauss_hermit(3, 2, RT);` 则是生成 D2Q9 ，且气体常数和温度乘积为 `RT` 的网格。

配置文件参数为：

```text
case-info:
    # other definitions
    ...
    # gauss-hermit
    dvs_mesh        GaussHermit
:end

def:
    # other definitions
    ...
    # gauss-hermit
    gauss_point     3
:end
```

调用代码如下：

```c++
ConfigReader config("./config/<your_config_file>");

/// other codes

// DVS mesh
MESH::StaticMesh dvs_mesh = MESH::StaticMesh(MESH_TYPE_NO_FACE, "GaussHermit");

D = phy_mesh.dimension();               // dimension
R = config.get<double>("R");            // gas constant
T = config.get<double>("temperature");  // temperature

if (config.dvs_mesh == MESH_GAUSS_HERMIT) {
        int gp;
        gp = config.get<int>("gauss_point");
        dvs_mesh = GENERATOR::gauss_hermit(gp, D, R * T);
    }
```

### Newton-Cotes型DVS

生成 Newton-Cotes 型结构速度空间网格的函数。其中 `n + 1` 和 `mount` 分别为复化求积中小区间内的结点个数和小区间份数。

如： `dvs_mesh = GENERATOR::newton_cotes(4, 22, 2, 15.0, RT);` 则是生成点数为 `89 * 89` ，范围为 `[-15*sqrt(2*RT), 15*sqrt(2*RT)] *  [-15*sqrt(2RT), 15*sqrt(2RT)]` 的二维结构速度空间，复化求积的代数精度为 7 阶。

配置文件参数为：

```text
case-info:
    # other definitions
    ...
    # newton-cotes
    dvs_mesh                NewtonCotes
:end

def:
    # other definitions
    ...
    # newton-cotes
    newton_cotes_point      4
    newton_cotes_mount      22
    newton_cotes_scale      15.0
:end
```

调用代码如下：

```c++
ConfigReader config("./config/<your_config_file>");

/// other codes

// DVS mesh
MESH::StaticMesh dvs_mesh = MESH::StaticMesh(MESH_TYPE_NO_FACE, "NewtonCotes");

D = phy_mesh.dimension();               // dimension
R = config.get<double>("R");            // gas constant
T = config.get<double>("temperature");  // temperature

if (config.dvs_mesh == MESH_NEWTON_COTES) {
        int n, mount;
        double scale;
        n = config.get<int>("newton_cotes_point");
        mount = config.get<int>("newton_cotes_mount");
        scale = config.get<double>("newton_cotes_scale");
        dvs_mesh = GENERATOR::newton_cotes(n, mount, D, scale, R * T);
    }
```
