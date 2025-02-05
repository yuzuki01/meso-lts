#ifndef CORE_LOGGER_H
#define CORE_LOGGER_H

namespace MESO {
    class Printer;
    class Logger;
}


class MESO::Printer {
private:
    Logger* root_ptr;
    const int output_level;
    const String color_prefix;
    const String color_suffix;
    bool isOutputAvailable();

public:
    explicit Printer(const int output_level, Logger* root_ptr, String prefix) :
            output_level(output_level),
            root_ptr(root_ptr),
            color_prefix(std::move(prefix)),
            color_suffix("\033[0m") {
    };

    template<class MesoType>
    Printer& operator<<(const MesoType& value) {
        if (isOutputAvailable()) {
            std::cout << color_prefix << value;
        }
        return *this;
    }

    Printer& operator<<(const Vector& value);

    Printer& operator<<(std::ostream& (*manipulator)(std::ostream&)) {
        if (isOutputAvailable()) {
            std::cout << color_suffix;
            manipulator(std::cout);
        }
        return *this;
    }
};


class MESO::Logger {
private:

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

namespace MESO {
    extern Logger logger;
}

#endif //CORE_LOGGER_H
