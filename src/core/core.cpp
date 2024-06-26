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
        logger.info << std::setw(15) << std::right << it.c_str();
    }
    logger.info << std::endl;
    for (auto &it: values) {
        logger.info << std::scientific << std::setw(15) << std::right << std::setprecision(4) << it;
    }
    logger.info << std::endl;
}

int Utils::mkdir(const MESO::String &dir_name) {
#ifdef _WIN32
    std::string command = "if not exist \"" + dir_name + "\" mkdir \"" + dir_name + "\"";
#else
    std::string command = "[ -d \"" + dir_name + "\" ] || mkdir \"" + dir_name + "\"";
#endif
    return system(command.c_str());
}

bool Utils::is_converged(const std::vector<double> &residual_list, double limit) {
    bool result = true;
    for (auto &it : residual_list) {
        result = result and (it <= limit);
    }
    return result;
}
