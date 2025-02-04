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
#include <vector>
#include <unordered_map>
#include <utility>
#include <random>
#include <memory>
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

    void print_names_and_values(const StringList &names, const List<Scalar> &values);

    int mkdir(const String &dir_name);

    bool is_converged(const std::vector<double> &residual_list, double limit);

    template<class T>
    void output_list(const String &file_path, std::vector<T> &data);

    template<class T>
    std::vector<T> read_np_file(const String &file_path);

    /// receive python3 script prints
    StringList exec_script(const String &script);
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
