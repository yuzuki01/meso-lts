#include "core/core.h"


namespace MESO {
    void init(MESO::ArgParser &parser);
}

using namespace MESO;

Logger logger;

StringList MESO::Utils::split(const std::string &_str) {
    std::istringstream iss(_str);
    std::vector<std::string> data;
    std::string token;
    while (iss >> token) {
        data.push_back(token);
    }
    return data;
}
