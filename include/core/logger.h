#ifndef CORE_LOGGER_H
#define CORE_LOGGER_H

namespace MESO {
    class Logger;
}


class MESO::Logger {
private:
    class Printer {
    private:
        Logger* root_ptr;
        const int output_level;
        const std::string color_prefix;
        const std::string color_suffix;
    public:
        explicit Printer(const int output_level, Logger* root_ptr, std::string prefix, std::string suffix="\033[0m") :
                output_level(output_level),
                root_ptr(root_ptr),
                color_prefix(std::move(prefix)),
                color_suffix(std::move(suffix)) {};

        template<typename T>
        Printer& operator<<(const T& value) {
            if (root_ptr->level <= output_level) {
                std::cout << color_prefix << value << color_suffix;
            }
            return *this;
        }

        Printer& operator<<(std::ostream& (*manipulator)(std::ostream&)) {
            if (root_ptr->level <= output_level) {
                manipulator(std::cout);
            }
            return *this;
        }
    };

public:
    /// indicate printer
    int level = 0;
    Printer info, note, warn, error, debug;

    Logger() :
            info(0, this, "\033[0m"),
            note(1, this, "\033[1;36m"),
            warn(2, this, "\033[1;31m"),
            error(3, this, "\033[41;37m"),
            debug(-1, this, "\033[1;35m") {};
};

extern MESO::Logger logger;

#endif //CORE_LOGGER_H
