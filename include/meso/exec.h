#ifndef MESO_MPI_EXEC_H
#define MESO_MPI_EXEC_H

#include <iostream>
#include <string>
#include <ctime>

// 预处理指令适配不同编译器，提取编译环境信息
#if defined(_MSC_VER)  // MSVC (Visual Studio)
#define COMPILER_NAME "MSVC"
#define COMPILER_VERSION _MSC_FULL_VER  // 完整版本号（如 193933520）
#define COMPILER_VERSION_STR _MSC_VER  // 主版本号（如 1930）
#elif defined(__GNUC__)  // GCC/G++
#define COMPILER_NAME "GCC/G++"
#define COMPILER_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#define COMPILER_VERSION_STR std::to_string(__GNUC__) + "." + std::to_string(__GNUC_MINOR__) + "." + std::to_string(__GNUC_PATCHLEVEL__)
#elif defined(__clang__)  // Clang/LLVM
#define COMPILER_NAME "Clang"
#define COMPILER_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#define COMPILER_VERSION_STR std::to_string(__clang_major__) + "." + std::to_string(__clang_minor__) + "." + std::to_string(__clang_patchlevel__)
#else
#define COMPILER_NAME "Unknown Compiler"
#define COMPILER_VERSION 0
#define COMPILER_VERSION_STR "Unknown"
#endif

// 编译时间宏处理（__DATE__ 和 __TIME__ 是编译器内置宏）
#define COMPILE_DATE __DATE__
#define COMPILE_TIME __TIME__


namespace MESO {
    void help();

    void version();

    int handle_mesh(const std::string &mesh_file, double mesh_scale = 1.0);
}

#endif //MESO_MPI_EXEC_H
