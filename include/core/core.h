#ifndef MESO_CORE_H
#define MESO_CORE_H

#include <ctime>
#include <cmath>
#include <csignal>
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

/// UDF
#include "core/math.h"
#include "core/type_def.h"
#include "core/mesompi.h"
#include "core/logger.h"
#include "core/argparser.h"

namespace MESO::DefaultValue {
    const int max_step = 1000000;
    const int save_interval = 1000;
}

namespace MESO::Utils {
    StringList split(const std::string &_str);

    void print_names_and_values(const StringList &names, const ScalarList &values);

    int mkdir(const String &dir_name);

    bool is_converged(const std::vector<double> &residual_list, double limit);
}

namespace MESO::FileReader {
    class BasicReader {
    public:
        String file;
        explicit BasicReader(String file_path) : file(std::move(file_path)) {
            logger.debug << "Read file: " << file << std::endl;
        };

        [[nodiscard]] StringList read_lines() const;
    };
}

#endif //MESO_CORE_H
