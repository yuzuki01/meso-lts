#ifndef CORE_ARGPARSER_H
#define CORE_ARGPARSER_H

#define PARSER_NULL "NULL"

namespace MESO {
    class ArgParser;
}

class MESO::ArgParser {
protected:
    std::unordered_map<std::string, bool> _switch;
    std::unordered_map<std::string, std::string> _value;    /// for param parser
public:
    ArgParser(int argc, char **argv);

    bool parse_switch(const std::string &target);

    template<class T>
    T parse_param(const std::string &target, const T &_default, bool print_error=false);
};

#endif //CORE_ARGPARSER_H
