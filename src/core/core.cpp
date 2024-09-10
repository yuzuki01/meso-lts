#include "core/core.h"


namespace MESO {
    void init(MESO::ArgParser &parser);
}

using namespace MESO;

Logger logger;

StringList Utils::split(const String &_str) {
    std::istringstream iss(_str);
    StringList data;
    String token;
    while (iss >> token) {
        token.erase(std::remove(token.begin(), token.end(), '\r'), token.end());
        token.erase(std::remove(token.begin(), token.end(), '\n'), token.end());
        data.push_back(token);
    }
    return data;
}

void Utils::print_names_and_values(const StringList &names, const ScalarList &values) {
    if (names.size() != values.size()) {
        throw std::invalid_argument("MESO::Utils::print_names_and_values() got unmatched size.");
    }
    for (auto &it: names) {
        logger.info << std::setw(10) << std::left << it.c_str() << "\t";
    }
    logger.info << std::endl;
    for (auto &it: values) {
        logger.info<< std::setw(10) << std::left << std::scientific << std::setprecision(2) << it << "\t";
    }
    logger.info << std::endl;
}

int Utils::mkdir(const MESO::String &dir_name) {
    if (MPI::rank != MPI::main_rank) return 0;
    std::string command = "[ -d \"" + dir_name + "\" ] || mkdir \"" + dir_name + "\"";
    return system(command.c_str());
}

bool Utils::is_converged(const std::vector<double> &residual_list, double limit) {
    bool result = true;
    for (auto &it : residual_list) {
        result = result and (it <= limit);
    }
    return result;
}
