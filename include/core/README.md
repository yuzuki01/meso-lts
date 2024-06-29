# Core

## Vector

3 维向量

| 运算 |                             重载运算符                             |
|:--:|:-------------------------------------------------------------:|
| 加法 |                            `a + b`                            |
| 减法 |                            `a - b`                            |
| 点乘 |                            `a * b`                            |
| 叉乘 |                            `a ^ b`                            |
| 比例 | `k * a` <br> `a * k` <br> `a / k` <br> `a *= k` <br> `a /= k` |
| 自加 |                           `a += b`                            |
| 自减 |                           `a -= b`                            |

## MESO MPI

为了支持 `MESO::Field` 的分布式并行，自定义了适用于 `MESO::Field` 的广播和归约函数

全局变量和函数位于命名空间 `MESO::MPI` 中

### MPI Rank

```c++
using namespace MESO;

logger.info << "Hello from node-" << MPI::rank << std::endl;

```

### Broadcast for Field

```c++
using namespace MESO;

Field<Scalar> f_scalar_global;
MPI::Bcast(f_scalar_global);

Field<Vector> f_vector_global;
MPI::Bcast(f_vector_global);

```

### Reduction for Field

```c++
using namespace MESO;

Field<Scalar> f_scalar_local, f_scalar_global;
MPI::AllReduce(f_scalar_local, f_scalar_global);

Field<Vector> f_vector_local, f_vector_global;
MPI::AllReduce(f_vector_local, f_vector_global);

```
