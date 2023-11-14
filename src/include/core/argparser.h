/**
 * included by core.h
 */

#define PARSER_NULL "NULL"

class ArgParser {
protected:
    std::unordered_map<std::string, bool> _switch;
    std::unordered_map<std::string, std::string> _value;    /// for param parser
public:
    ArgParser(int argc, char **argv);

    bool parse_switch(const std::string &target);

    template<class T>
    T parse_param(const std::string &target, const T &_default);
};
