/**************************
 * Logger for meso        *
 *                        *
 **************************/

# include <core.h>

std::string logo = "******************************\n"
                   "*     __  ___                *\n"
                   "*    /  |/  /__  _________   *\n"
                   "*   / /|_/ / _ \\/ ___/ __ \\  *\n"
                   "*  / /  / /  __(__  ) /_/ /  *\n"
                   "* /_/  /_/\\___/____/\\____/   *\n"
                   "*                            *\n"
                   "******************************\n";

/**
 * logger 基类
 */

Logger::Logger(const std::string &logger_name) {
    ss.str("");
    path.clear();
    prefix = logger_name;
}

Logger::Logger(const std::string &logger_name, const std::string &file_path) {
    ss.str("");
    prefix = logger_name;
    path = file_path;
}

Logger &Logger::operator()(const std::string &logger_name) {
    ss.str("");
    path.clear();
    prefix = logger_name;
    return *this;
}

Logger &Logger::operator<<(const std::string &_text) {
    ss << _text;
    return *this;
}

Logger &Logger::operator<<(int _value) {
    ss << _value;
    return *this;
}

Logger &Logger::operator<<(double _value) {
    ss << _value;
    return *this;
}

void Logger::info() {
    std::cout << ColorMark::none << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::note() {
    std::cout << ColorMark::note << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::warn() {
    std::cout << ColorMark::warn << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::error() {
    std::cout << ColorMark::error << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::debug() {
    std::cout << ColorMark::debug << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::green() {
    std::cout << ColorMark::green << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}

void Logger::highlight() {
    std::cout << ColorMark::highlight << "[" << prefix << "] " << ss.str() << ColorMark::none << std::endl;
    ss.str("");
}


void LoggerInit() {
#if OS == 0
    system("cls");
#elif OS == 1
    system("clear");
#endif
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);
    std::cout << logo;
}
