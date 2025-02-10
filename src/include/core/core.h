#ifndef MESO_CORE_H
#define MESO_CORE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <csignal>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <utility>
#include <random>
#include <memory>
#include <mpi.h>


/// Global
namespace MESO {
    extern bool debug;

    extern const std::string logo;

    const int DATA_PRECISION = 18;
    const int LINE_DATA_NUM = 30;
}

/// MACRO
#define namespaceMESO using namespace MESO;

#define FATAL_ERROR_THROW  throw std::invalid_argument("MESO FATAL ERROR")

#define forAll(X, i) for(int i=0; i < (X).size(); ++i)

#define forConstRef(X, List) for(const auto& X : List)

/// Costume
#include "core/typeDefine.h"
#include "core/math.h"
#include "core/mesompi.h"
#include "core/logger.h"
#include "core/argparser.h"
#include "core/utils.h"
#include "core/fileReader.h"
#include "core/timeControl.h"

#endif //MESO_CORE_H
