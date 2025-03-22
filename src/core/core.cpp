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

void Utils::print_names_and_values(const StringList &names, const List<Scalar> &values) {
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

int Utils::mkdir(const String &dir_name) {
    if (MPI::rank != MPI::main_rank) return 0;
    std::stringstream ss;
    ss << "mkdir -p " << dir_name;
    return system(ss.str().c_str());
}

bool Utils::is_converged(const List<Scalar> &residual_list, Scalar limit) {
    bool result = true;
    for (auto &it : residual_list) {
        result = result and (it <= limit);
    }
    return result;
}

template<>
void Utils::output_list(const String &file_path, List<Scalar> &data) {
    if (MPI::rank != MPI::main_rank) return;
    const int DATA_PRECISION=18;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &it : data) {
        fp << std::setprecision(DATA_PRECISION) << it << std::endl;
    }
}

template<>
void Utils::output_list(const String &file_path, List<Vector> &data) {
    if (MPI::rank != MPI::main_rank) return;
    const int DATA_PRECISION=18;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &it : data) {
        fp << std::setprecision(DATA_PRECISION) << it.x << " " << it.y << " " << it.z << std::endl;
    }
}

template<>
void Utils::output_list(const String &file_path, List<List<ObjectId>> &data) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &list : data) {
        fp << "(";
        for (auto &it: list) {
            fp << " " << it;
        }
        fp << " )" << std::endl;
    }
}


template<>
void Utils::output_list(const String &file_path, List<Set<ObjectId>> &data) {
    if (MPI::rank != MPI::main_rank) return;
    std::fstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    for (auto &set : data) {
        fp << "( " << set[0] << " " << set[1] << " )" << std::endl;
    }
}

template<>
List<ObjectType> Utils::read_np_file<ObjectType>(const String &file_path) {
    List<ObjectType> result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.push_back(stoi(data[0]));
    }
    return result;
}

template<>
List<Scalar> Utils::read_np_file<Scalar>(const String &file_path) {
    List<Scalar> result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.push_back(stod(data[0]));
    }
    return result;
}

template<>
List<Vector> Utils::read_np_file<Vector>(const String &file_path) {
    List<Vector> result;
    MESO::FileReader::BasicReader reader(file_path);
    auto lines = reader.read_lines();
    for (auto &line: lines) {
        auto data = Utils::split(line);
        result.emplace_back(stod(data[0]), stod(data[1]), stod(data[2]));
    }
    return result;
}

StringList Utils::exec_script(const String &script) {
    std::array<char, 1024> buffer{};
    StringList result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
            popen(("python3 " + script).c_str(), "r"), pclose);

    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        std::string line(buffer.data());
        if (!line.empty() && line.back() == '\n') {
            line.erase(line.length() - 1);
        }
        if (!line.empty() && line.back() == '\r') {
            line.erase(line.length() - 1);
        }
        result.push_back(line);
    }

    return result;
}
