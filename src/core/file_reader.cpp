#include "core/core.h"


using namespace MESO::FileReader;


MESO::StringList BasicReader::read_lines() const {
    std::ifstream fp(file);
    if (!fp.is_open()) {
        logger.warn << "Cannot open the file: " << file << std::endl;
        fp.close();
        throw std::invalid_argument("file open error");
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
