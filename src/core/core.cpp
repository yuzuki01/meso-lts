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
    std::stringstream ss;
    ss << "mkdir -p " << dir_name;
    return system(ss.str().c_str());
}

bool Utils::is_converged(const std::vector<double> &residual_list, double limit) {
    bool result = true;
    for (auto &it : residual_list) {
        result = result and (it <= limit);
    }
    return result;
}

template<>
void Utils::output_list(const MESO::String &file_path, std::vector<Scalar> &data) {
    const int DATA_PRECISION=18;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &it : data) {
        fp << std::setprecision(DATA_PRECISION) << it << std::endl;
    }
}

template<>
void Utils::output_list(const MESO::String &file_path, std::vector<Vector> &data) {
    const int DATA_PRECISION=18;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &it : data) {
        fp << std::setprecision(DATA_PRECISION) << it.x << " " << it.y << " " << it.z << std::endl;
    }
}

template<>
ObjectTypeList Utils::read_np_file<ObjectType>(const MESO::String &file_path) {
    ObjectTypeList result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.push_back(stoi(data[0]));
    }
    return result;
}

template<>
ScalarList Utils::read_np_file<Scalar>(const MESO::String &file_path) {
    ScalarList result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.push_back(stod(data[0]));
    }
    return result;
}

template<>
VectorList Utils::read_np_file<Vector>(const MESO::String &file_path) {
    VectorList result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.emplace_back(stod(data[0]), stod(data[1]), stod(data[2]));
    }
    return result;
}
