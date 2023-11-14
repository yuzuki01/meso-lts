/**
 * included by core.h
 */

/// Global Vars
namespace ColorMark {
    using cst_string = const std::string;
    cst_string none = "\033[0m";
    cst_string note = "\033[1;36m";
    cst_string warn = "\033[1;31m";
    cst_string error = "\033[41;37m";
    cst_string debug = "\033[1;35m";
    cst_string green = "\033[1;32m";
    cst_string highlight = "\033[1;33m";
}

/// 初始化函数
void LoggerInit();

/// Logger 类
class Logger {
protected:
    std::string prefix{};
    std::string path{};
    std::stringstream ss{};
public:

    Logger() = default;

    Logger &operator()(const std::string &logger_name);

    explicit Logger(const std::string &logger_name);

    Logger(const std::string &logger_name, const std::string &file_path);

    Logger &operator<<(const std::string &_text);

    Logger &operator<<(int _value);

    Logger &operator<<(double _value);

    void info();

    void note();

    void warn();

    void error();

    void debug();

    void green();

    void highlight();
};
