/**************************************
 * meso  - Core module                *
 *  contains:                         *
 *      basic class, global functions *
 *                Nov 6, 2023  by MYC *
 **************************************/

#ifndef HEADER_CORE
#define HEADER_CORE

/**
 * 基础库
 */
#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <exception>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <omp.h>

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

/**
 * Global Vars
 */

extern bool debug_mode;

/**
 * Macro define
 */

#include "core/macro.h"

/**
 * Basic class
 */

#include "core/argparser.h"
#include "core/logger.h"
#include "core/mesomath.h"
#include "core/reader.h"

/**
 * Global func
 */

#include "core/core_utils.h"

#endif  // HEADER_CORE
