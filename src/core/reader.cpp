/**************************
 * Reader for meso        *
 *                        *
 **************************/

#include <core.h>

BasicReader::BasicReader(const std::string &reader_name, const std::string &file_path) {
    name = reader_name;
    path = file_path;
    logger(name);
    std::ifstream file(file_path);
    file_open = file.is_open();
    if (!file_open) {
        logger << "Warn: Cannot open the file " << file_path;
        logger.warn();
        file.close();
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        /// skip empty line
        if (!line.empty()) {
            /// skip #
            if (strncmp(line.c_str(), "# ", 2) != 0) {
                lines.push_back(line);
            }
        }
    }
    file.close();
    logger << "Open file: " << file_path;
    logger.note();
}

bool BasicReader::is_file_open() const {
    return file_open;
}

void BasicReader::print() {
    if (file_open) {
        logger << "Context of reader[" << name << "]:";
        logger.note();
        /// print context
        int i = 1;
        for (auto &line : lines) {
            std::stringstream ss;
            ss << "l" << std::setw(2) << i++;
            info_println(ss.str() + "| " + line);
            if (i > print_line_num && lines.size() > print_line_num) {
                ss.str("");
                ss << "...left " << lines.size() + 1 - i << " lines.";
                info_println(ss.str());
                break;
            }
        }
    } else {
        logger << "Warn: Opening file failed.";
        logger.warn();
    }
}

int BasicReader::line_num() const {
    return lines.size();
}

std::string & BasicReader::operator[](int index) {
    return lines[index];
}
