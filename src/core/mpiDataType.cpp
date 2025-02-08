#include "core/core.h"


namespaceMESO


/// Indication
MPI_Datatype MPI::DataType::MPI_VECTOR;


void MPI::DataType::MPI_Vector_init() {
    // 定义结构体中各成员的数据块长度
    int blocklengths[3] = {1, 1, 1};
    // 定义每个成员的 MPI 类型
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[3];

    // 使用一个 dummy 实例来计算偏移量
    Vector dummy;
    MPI_Aint base_address;
    MPI_Get_address(&dummy, &base_address);
    MPI_Get_address(&dummy.x, &displacements[0]);
    MPI_Get_address(&dummy.y, &displacements[1]);
    MPI_Get_address(&dummy.z, &displacements[2]);

    // 转换为相对于结构体起始地址的偏移量
    displacements[0] -= base_address;
    displacements[1] -= base_address;
    displacements[2] -= base_address;

    // 创建并提交 MPI 自定义数据类型
    MPI_Type_create_struct(3, blocklengths,
                           displacements, types,
                           &DataType::MPI_VECTOR);
    MPI_Type_commit(&DataType::MPI_VECTOR);
}

void MPI::DataType::MPI_Vector_free() {
    MPI_Type_free(&DataType::MPI_VECTOR);
}
