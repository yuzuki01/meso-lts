# 求解器配置文件

## settings

```
[settings]
case-name   <your-project-name>
mesh-file   <your-mesh-file>
dvs-file    <your-dvs-file>
dvs-type    <your-dvs-type>
```

### 必备参数

`case-name` 会在可执行文件所在目录下创建同名文件夹用于储存结果

`mesh-file` 可以填写绝对路径或相对路径

`dvs-file`  非结构速度空间网格路径或生成结构网格的配置文件路径

### 可选参数

`mesh-scale` 用于放缩网格尺寸，默认 `1.0`，不进行放缩

`dvs-type`  决定速度空间网格类型，生成网格的脚本位于 `files/scritps`

 - `Gauss-Hermite`
 - `Newton-Cotes`
 - `<dvs-type>`：默认值，或其他类型，都调用 `fvmMesh::load_gambit()` 函数

继承于基类 `MESO::Solver::BasicSolver` 的参数：

`case-name`：算例名称，默认 `unnamed`

`output-np`：是否打开 `numpy` 格式字符输出，默认`False`

`read-np-data`：是否读取 `numpy` 格式输出作为初场，默认 `False`

`residual-limit`：残差限，默认 `1e-6`，在所有监测量的残差连续 `10` 次监测都小于残差限后，即使未达到最大迭代步数求解器也会停止并保存。

---

根据求解器内编程定义（具体定义根据求解器代码而定）：

对于 `cdugks@incompressible` ，已经定义有:

流体力学参数：

`Re`，`Ma`，`CFL`

求解器设置：

`limiter-switch`：是否打开限制器，默认`False`

 - `venkata-limiter-k`：限制器参数，默认`1.0`，设置为`0`时限制效果最强

`gradient-switch`：是否打开梯度计算，默认`True`

对于 `cdugks@shakhov` ，已经定义有：

`dvs-file`：离散速度空间网格文件，填写 `Newton-Cotes` 则采用牛顿科特斯网格，自动生成。

## group

网格内划分的单元区域

```
[group]
name        <group-name>
density     <value>
temperature <value>
veloxity-x  <value>
...
```

## 边界定义

|      type       |  边界含义  |                   物理量                   |
|:---------------:|:------:|:---------------------------------------:|
| farfield-inlet  |  远场入口  |    **给定密度**<br>**给定温度**<br>**给定速度**     |
| pressure-inlet  |  压力入口  |    **给定密度**<br>**给定温度**<br>_邻近单元速度_     |
| farfield-outlet |  远场出口  |                _邻近单元赋值_                 |
| pressure-outlet |  压力出口  |    **给定密度**<br>**给定温度**<br>_邻近单元速度_     |
|      wall       |  等温固壁  | **给定温度**<br>**计算无穿透密度**<br>_Maxwell 扩散_ |
|    slip-wall    | 等温滑移固壁 |            Maxwell 扩散 + 镜面反射            |
|    symmetry     |  对称边界  |               分布函数按界面对称赋值               |
|      patch      | 自定义边界  |     通过定义 `PatchType` 实现不同的边界变量计算方式      |

### Patch

**PatchType**

`calculated`: 标记为需要计算，具体算法在代码中实现

`fixedValue`: 标记为定值，类型包含 `scalar` 和 `vector`

* 示例: `density    fixedValue    scalar    1.0`

`zeroGradient`: 标记为临近单元赋值，对于界面使用其临近单元的值进行赋值

`fromFile`: 标记为从 `.np.dat` 格式文件读取，按照网格对象的编号进行索引，类型包含 `scalar` 和 `vector`

* 示例: `velocity    fromFile    vector    vel.np.dat`

```
[mark]
name            <name>
type            <mark-type>
density         fixedValue    scalar    1.0
velocity        fixedValue    vector    1.0 0.0 0.0
temperature     fromFile      scalar    inlet-T.np.dat
...
```
