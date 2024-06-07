# 求解器配置文件

## settings

```
[settings]
case-name   <your-project-name>
mesh-file   <your-mesh-file>
```

### 必备参数

`case-name` 会在可执行文件所在目录下创建同名文件夹用于储存结果

`mesh-file` 可以填写绝对路径或相对路径

### 可选参数

根据求解器内编程定义

对于 `不可压缩 cdugks` ，已经定义有:

流体力学参数：

`Re`，`Ma`，`CFL`

求解器设置：

`limiter-switch`：是否打开限制器，默认`False`

 - `venkata-limiter-k`：限制器参数，默认`1.0`，设置为`0`时限制效果最强

`gradient-switch`：是否打开梯度计算，默认`True`

`output-np`：是否打开 `numpy` 格式字符输出，默认`False`

## group

网格内划分的单元区域

```
[group]
name        <group-name>
density     <value>
veloxity-x  <value>
...
```

边界定义

```
[mark]
name        <name>
type        <mark-type>
density     <value>
```
