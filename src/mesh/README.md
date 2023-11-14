# meso/mesh

meso 网格

## 模块

### 对象定义

当前网格分类固定网格（MESH::ListMesh，容器 std::vector）和可变网格（MESH::MapMesh 容器 std::unordered_map）

可变网格支持网格几何对象的增减

不同的网格类通过模板函数实现了通用接口

```c++
#define ListMesh Mesh<int>
#define MapMesh Mesh<std::string>

MESH::ListMesh list_mesh;
MESH::MapMesh map_mesh;
```

其对应的网格几何对象也不同

```c++
#define TP_key template <class key_type>
#define std::vector<key_type> key_vector
#define std::vector<std::string> string_vector
/// 网格结点
TP_key
class MESH::Node {
public:
    // key 值
    key_type key;
    // 坐标
    Vec3D position = {0.0, 0.0, 0.0};
    // 按照 cpp 结构构造
    explicit Node(const key_type &node_key, const Vec3D &node_position);
    // 按照 Gambit 格式构造
    explicit Node(const string_vector &init_strings);
};

/// 网格单元
TP_key
class MESH::Cell {
public:
    // Gambit 几何体类型
    int type;
    // key 值
    key_type key;
    // 对于结构速度空间 - 反方向 cell 的 key 值
    key_type inv_cell_key;
    // 对于结构速度空间 - 镜面反射
    int on_layer, layer_num;
    // 几何参数
    // 界面个数
    int face_num;
    // 影单元
    bool have_shadow;
    key_type shadow_cell_key;
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
    /// 按照 Gambit 格式构造
    explicit Cell(const string_vector &init_strings);
    /// 结构速度空间构造函数
    Cell(const Vec3D &particle_velocity, double weight,
         const key_type &cell_key, const key_type &inv_cell_key,
         int on_layer_cell_id, int on_layer_cell_num);
};

/// 网格界面
TP_key
class MESH::Face {
public:
    // Gambit 几何体类型
    int type;
    // key 值
    key_type key;
    // 边界类型     默认为 0 ，非边界界面
    int boundary_type = 0;
    // 界面所在边界(MESH::Mark) key 值
    std::string mark_key = MESH_KEY_NULL;

    /// 几何参数
    /**
     * on_cell 表示正向 inv_cell 表示反向
     *  界面为边界时，由于只有一个相邻单元(MESH::Cell)
     *  此时：
     *      on_cell_key == inv_cell_key
     *      on_cell_face == inv_cell_face
     *      on_cell_nv == -inv_cell_nv
     *      on_cell_nv 为单元格心向外的界面法向量
     *      inv_cell_nv 为边界指向流体域的界面法向量
     */
    key_type on_cell_key;
    int on_cell_face;
    Vec3D on_cell_nv;
    key_type inv_cell_key;
    int inv_cell_face;
    Vec3D inv_cell_nv;
    // 面积
    double area;
    // 坐标
    Vec3D position = {0.0, 0.0, 0.0};
    // 构成界面的节点(MESH::Node) key 值
    key_vector node_key;

    /// 构造
    Face(const key_type &face_key, int elem_type, const key_vector &node_key_vec);
    Face();
    
    // get item
    Node<key_type> &get_node(const key_type &_key);

    Cell<key_type> &get_cell(const key_type &_key);

    Face<key_type> &get_face(const key_type &_key);

    Mark<key_type> &get_mark(const key_type &_key);
};
```

## 网格生成器

### 生成 Gauss-Hermit 型 DVS
```c++
ConfigReader config("./config/<your_config_file>");

/// other codes

// DVS mesh
MESH::Mesh<int> dvs_mesh = MESH::MapMesh(MESH_TYPE_NO_FACE);

D = phy_mesh.dimension();               // dimension
R = config.get<double>("R");            // gas constant
T = config.get<double>("temperature");  // temperature

if (config.dvs_mesh == MESH_GAUSS_HERMIT) {
        int gp;
        gp = config.get<int>("gauss_point");
        dvs_mesh = GENERATOR::gauss_hermit(gp, D, R * T);
    }
```
