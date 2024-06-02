#ifndef MESO_CORE_H
#define MESO_CORE_H

#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <exception>
#include <array>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <random>
#include <functional>
#include <omp.h>
#include <mpi.h>

#if defined(_WIN32)
/// Windows
#define OS 0
#include <direct.h>
#elif defined(__linux__)
/// Linux
#define OS 1
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#else
/// Unsupported platform
#define OS -1
#endif


/// UDF
#include "core/math.h"
#include "core/type_def.h"
#include "core/mesompi.h"
#include "core/logger.h"
#include "core/argparser.h"

namespace MESO::Utils {
    StringList split(const std::string &_str);
}

#endif //MESO_CORE_H
