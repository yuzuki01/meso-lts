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
        bool is_output_available() {
            return (root_ptr->level <= output_level) and (MPI::rank == MPI::main_rank);
        };

    public:
        explicit Printer(const int output_level, Logger* root_ptr, std::string prefix) :
                output_level(output_level),
                root_ptr(root_ptr),
                color_prefix(std::move(prefix)) {};

        template<typename T>
        Printer& operator<<(const T& value) {
            if (is_output_available()) {
                std::cout << color_prefix << value;
            }
            return *this;
        }

        Printer& operator<<(std::ostream& (*manipulator)(std::ostream&)) {
            if (is_output_available()) {
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
