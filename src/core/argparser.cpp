/**************************
 * ArgParser              *
 *      for meso          *
 *                        *
 **************************/

#include "core.h"

ArgParser::ArgParser(int argc, char **argv) {
    for (int i = 0; i < argc; i++) {
        if (std::strncmp(argv[i], "-", 1) == 0) {
            if (std::strncmp(argv[i], "--", 2) == 0) {
                /// --target value
                _value[argv[i]] = (argc > i + 1 && std::strncmp(argv[i + 1], "-", 1) != 0) ? argv[i + 1] : PARSER_NULL;
            } else {
                /// -target
                _switch[argv[i]] = true;
            }
        }
    }
}

bool ArgParser::parse_switch(const std::string &target) {
    auto it = _switch.find("-" + target);
    if (it == _switch.end()) {
        /// not found
        return false;
    } else {
        /// found
        return _switch["-" + target];
    }
}


template<> std::string ArgParser::parse_param<std::string>(const std::string &target,
        const std::string &_default) {
    auto it = _value.find("--" + target);
    if (it == _value.end()) {
        return _default;
    } else {
        return _value["--" + target];
    }
}

template<> bool ArgParser::parse_param<bool>(const std::string &target, const bool &_default) {
    auto it = _value.find("--" + target);
    if (it == _value.end()) return _default;
    auto &it_value = _value["--" + target];
    if (it_value == "true" || it_value == "True") return true;
    return false;
}

template<> int ArgParser::parse_param<int>(const std::string &target,
        const int &_default) {
    auto it = _value.find("--" + target);
    if (it == _value.end()) {
        return _default;
    } else {
        return stoi(_value["--" + target]);
    }
}

template<> double ArgParser::parse_param<double>(const std::string &target,
                                              const double &_default) {
    auto it = _value.find("--" + target);
    if (it == _value.end()) {
        return _default;
    } else {
        return stod(_value["--" + target]);
    }
}
