#include "core/core.h"

using namespace MESO;
using namespace FileIO;

/**
 * =======================================================
 * ------------------- Basic Reader ----------------------
 * =======================================================
 **/

BasicReader::BasicReader(const String &filePath) : file_(filePath) {
    logger.debug << "Read file: " << filePath << std::endl;
}


StringList BasicReader::read_lines() const {
    std::ifstream fp(file_);
    if (!fp.is_open()) {
        logger.error << "Cannot open the file: " << file_ << std::endl;
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
    return file_;
}


/**
 * =======================================================
 * ------------------- Param Reader ----------------------
 * =======================================================
 **/

ParamReader::ParamReader(const String &filePath)
        : BasicReader(filePath) {
    update();
}

bool ParamReader::isVarExisted(const String &varRegion, const String &varName) {
    if (data_.find(varRegion) == data_.end()) return false;
    if (data_[varRegion].find(varName) == data_[varRegion].end()) return false;
    return true;
}

void ParamReader::clear() {
    data_.clear();
}

template<>
void ParamReader::set(const String &varRegion, const String &varName, const String &varValue) {
    data_[varRegion][varName] = varValue;
}

template<>
void ParamReader::set(const String &varRegion, const String &varName, const Label &varValue) {
    data_[varRegion][varName] = std::to_string(varValue);
}

template<>
void ParamReader::set(const String &varRegion, const String &varName, const Scalar &varValue) {
    data_[varRegion][varName] = std::to_string(varValue);
}

template<>
void ParamReader::set(const String &varRegion, const String &varName, const Vector &varValue) {
    StringStream ss;
    ss << varValue.x << " " << varValue.y << " " << varValue.z;
    data_[varRegion][varName] = ss.str();
}

template<>
String ParamReader::get(const String &varRegion, const String &varName, const String &varDefault,
                        bool throwNotFoundErr) {
    if (isVarExisted(varRegion, varName)) {
        return data_[varRegion][varName];
    } else if (throwNotFoundErr) {
        logger.warn << "ParamReader[" << name() << "] tried to fetch a nonexistent value: "
                    << varRegion << "::" << varName << std::endl;
        FATAL_ERROR_THROW;
    }
    return varDefault;
}

template<>
Label ParamReader::get(const String &varRegion, const String &varName, const Label &varDefault,
                       bool throwNotFoundErr) {
    if (isVarExisted(varRegion, varName)) {
        return std::stoi(data_[varRegion][varName]);
    } else if (throwNotFoundErr) {
        logger.warn << "ParamReader[" << name() << "] tried to fetch a nonexistent value: "
                    << varRegion << "::" << varName << std::endl;
        FATAL_ERROR_THROW;
    }
    return varDefault;
}


template<>
Scalar ParamReader::get(const String &varRegion, const String &varName, const Scalar &varDefault,
                        bool throwNotFoundErr) {
    if (isVarExisted(varRegion, varName)) {
        return std::stod(data_[varRegion][varName]);
    } else if (throwNotFoundErr) {
        logger.warn << "ParamReader[" << name() << "] tried to fetch a nonexistent value: "
                    << varRegion << "::" << varName << std::endl;
        FATAL_ERROR_THROW;
    }
    return varDefault;
}


template<>
Vector ParamReader::get(const String &varRegion, const String &varName, const Vector &varDefault,
                        bool throwNotFoundErr) {
    if (isVarExisted(varRegion, varName)) {
        auto split = Utils::split(data_[varRegion][varName]);
        return {std::stod(split[0]), std::stod(split[1]), std::stod(split[2])};
    } else if (throwNotFoundErr) {
        logger.warn << "ParamReader[" << name() << "] tried to fetch a nonexistent value: "
                    << varRegion << "::" << varName << std::endl;
        FATAL_ERROR_THROW;
    }
    return varDefault;
}

void ParamReader::update() {
    clear();
    const auto &lines = read_lines();
    String varRegion = "0";
    auto extractContent = [](const String &input) {
        size_t start_pos = input.find('[');
        size_t end_pos = input.find(']');
        if (start_pos != std::string::npos && end_pos != std::string::npos && start_pos < end_pos) {
            return input.substr(start_pos + 1, end_pos - start_pos - 1);
        } else {
            return fileAnnotation;
        }
    };
    for (const auto &line: lines) {
        auto data = Utils::split(line);
        if (data[0] == fileAnnotation) continue;
        {
            auto tmp = extractContent(data[0]);
            if (tmp != fileAnnotation) {
                if (tmp.empty()) {
                    logger.error << "ParamReader caught unexpected format" << std::endl;
                    FATAL_ERROR_THROW;
                }
                varRegion = tmp;
                continue;
            }
        }
        if (data.size() > 1) {
            StringStream var;
            var << data[1];
            if (data.size() > 2)
                for (int i = 2; i < data.size(); ++i) {
                    var << " " << data[i];
                }
            data_[varRegion][data[0]] = var.str();
        }
    }
}
