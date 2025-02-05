#include "core/core.h"

using namespace MESO;
using namespace MESO::FileReader;


StringList BasicReader::read_lines() const {
    std::ifstream fp(file_path);
    if (!fp.is_open()) {
        logger.error << "Cannot open the file: " << file_path << std::endl;
        fp.close();
        FATAL_ERROR_THROW;
    }
    String line;
    StringList lines;
    while (std::getline(fp, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        /// skip empty line
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    fp.close();
    return lines;
}

const String &BasicReader::name() const {
    return file_path;
}
