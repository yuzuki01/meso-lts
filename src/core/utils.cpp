#include "core.h"

/**
 * Global Vars
 */

bool debug_mode = false;

/**
 * Global Funcs
 */

void meso_init() {
    /// Logger
    LoggerInit();
    /// File
    std::vector<std::string> dirs = {"config", "mesh", "result"};
    for (auto &it : dirs) {
        if (!create_dir(it)) {
            std::stringstream ss;
            ss << "Cannot create dir: " << it;
            warn_println(ss.str());
        }
    }
}

string_vector split(const std::string &_str) {
    std::istringstream iss(_str);
    std::vector<std::string> data;
    std::string token;
    while (iss >> token) {
        data.push_back(token);
    }
    return data;
}

bool str_vec_cmp(const string_vector &_str_vec, const std::string &reference) {
    string_vector _ref = split(reference);
    int i = _str_vec.size(), j = _ref.size();
    if (i < j) return false;
    for (i = 0; i < j; i++) {
        if (_str_vec[i] != _ref[i]) return false;
    }
    return true;
}


/**
 * print
 **/

void info_println(const std::string &context) {
    std::cout << ColorMark::none << context << ColorMark::none << std::endl;
}

void note_println(const std::string &context) {
    std::cout << ColorMark::note << context << ColorMark::none << std::endl;
}

void warn_println(const std::string &context) {
    std::cout << ColorMark::warn << context << ColorMark::none << std::endl;
}

void error_println(const std::string &context) {
    std::cout << ColorMark::error << context << ColorMark::none << std::endl;
}

void debug_println(const std::string &context) {
    std::cout << ColorMark::debug << context << ColorMark::none << std::endl;
}

void green_println(const std::string &context) {
    std::cout << ColorMark::green << context << ColorMark::none << std::endl;
}

void highlight_println(const std::string &context) {
    std::cout << ColorMark::highlight << context << ColorMark::none << std::endl;
}

void data_int_println(const std::vector<std::string> &names, const std::vector<int> &values, const int width) {
    if (names.size() != values.size()) {
        std::string s = "Caught lists with different sizes to output.";
        warn_println(s);
        return;
    }
    std::stringstream ssn, ssv;
    ssn << std::right;
    for (auto &var : names) ssn << std::setw(width) << var;
    ssv << std::right;
    for (auto &var : values) ssv << std::scientific << std::setprecision(OUT_PRECISION) << std::setw(width) << var;
    info_println(ssn.str());
    info_println(ssv.str());
}

void data_double_println(const std::vector<std::string> &names, const std::vector<double> &values, const int width) {
    if (names.size() != values.size()) {
        std::string s = "Caught lists with different sizes to output.";
        warn_println(s);
        return;
    }
    std::stringstream ssn, ssv;
    ssn << std::right;
    for (auto &var : names) ssn << std::setw(width) << var;
    ssv << std::right;
    for (auto &var : values) ssv << std::setw(width) << var;
    info_println(ssn.str());
    info_println(ssv.str());
}

void data_sci_double_println(const std::vector<std::string> &names, const std::vector<double> &values, const int width) {
    if (names.size() != values.size()) {
        std::string s = "Caught lists with different sizes to output.";
        warn_println(s);
        return;
    }
    std::stringstream ssn, ssv;
    ssn << std::right;
    for (auto &var : names) ssn << std::setw(width) << var;
    ssv << std::right;
    for (auto &var : values) ssv << std::setw(width) << std::scientific << var;
    info_println(ssn.str());
    info_println(ssv.str());
}

/**
 * file
 */
bool create_dir(const std::string &path) {
#if OS == 1   /// linux
    if (access(path.c_str(), 0) == -1) mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    return (access(path.c_str(), 6) == 0);
#elif OS == 0 /// windows
    if (_access(path.c_str(), 0) == -1) _mkdir(path.c_str());
    return (_access(path.c_str(), 6) == 0);
#else
    return false;
#endif
}
