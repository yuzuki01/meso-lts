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
    forAll(lines, li) {
        auto line = lines[li];
        auto data = Utils::split(line);
        if (data[0] == fileAnnotation) continue;
        {
            auto tmp = extractContent(data[0]);
            if (tmp != fileAnnotation) {
                if (tmp.empty()) {
                    logger.error << "ParamReader caught unexpected format" << std::endl;
                    FATAL_ERROR_THROW;
                }
                if (tmp == "mark" or tmp == "zone") {
                    li++;
                    PatchParam patch;
                    while (li < lines.size()) {
                        line = lines[li];
                        data = Utils::split(line);
                        if (data[0] == fileAnnotation) {
                            li++;
                            continue;
                        } else if (extractContent(data[0]) != fileAnnotation) {
                            li--;
                            break;
                        }
                        patch.set(line, li);
                        li++;
                    }
                    if (tmp == "mark") marks_[patch.name()] = patch;
                    if (tmp == "zone") zones_[patch.name()] = patch;
                    continue;
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
    marks_["fluidInterior"] = PatchParam("fluidInterior", "fluidInterior");
    marks_["processor"] = PatchParam("processor", "processor");
}

const ParamReader::PatchParam &ParamReader::zone(const String &name){
    return zones_[name];
}

const ParamReader::PatchParam &ParamReader::mark(const String &name){
    return marks_[name];
}

/**
 * =======================================================
 * ------------------- Param Reader ----------------------
 * -------------------   PatchParm  ----------------------
 * =======================================================
 **/

ParamReader::PatchParam::PatchParam(const String &name, const String &type) {
    dataStr_["name"] = name;
    dataStr_["type"] = type;
}

const String &ParamReader::PatchParam::name() const {
    return dataStr_.at("name");
}

const String &ParamReader::PatchParam::type() const {
    return dataStr_.at("type");
}

void ParamReader::PatchParam::set(const String &line, Label lineNo) {
    auto data = Utils::split(line);
    if (data[0] == "name") {
        dataStr_["name"] = data[1];
    } else if (data[0] == "type") {
        dataStr_["type"] = data[1];
    } else {
        if (data[0] == "#") return;
        if (data.size() < 3) {
            logger.warn << "ParamReader::PatchParam caught invalid format data:" << std::endl << "\t" << lineNo + 1 << "| ";
            forConstRef(it, data) {
                logger.info << it << " ";
            }
            logger.info << std::endl;
            FATAL_ERROR_THROW;
        }
        auto varName = data[0];
        auto varType = data[1];
        if (Flag::TypeName.find(varType) == Flag::TypeName.end()) {
            logger.warn << "ParamReader::PatchParam caught invalid format data:" << std::endl << "\t" << lineNo + 1 << "| ";
            forConstRef(it, data) {
                logger.info << it << " ";
            }
            logger.info << std::endl;
            FATAL_ERROR_THROW;
        }
        switch (Flag::TypeName.at(varType)) {
            case Flag::string:
                dataStr_[varName] = data[2];
                break;
            case Flag::scalar:
                dataScl_[varName] = std::stod(data[2]);
                break;
            case Flag::vector:
                if (data.size() < 5) {
                    logger.warn << "ParamReader::PatchParam caught invalid format data:" << std::endl << "\t" << lineNo + 1 << "| ";
                    forConstRef(it, data) {
                        logger.info << it << " ";
                    }
                    logger.info << std::endl;
                    FATAL_ERROR_THROW;
                }
                dataVtr_[varName] = Vector(
                        std::stod(data[2]),
                        std::stod(data[3]),
                        std::stod(data[4])
                );
                break;
            default:
                logger.warn << "PatchParam caught invalid data:" << std::endl << "\t";
                forConstRef(it, data) {
                    logger.warn << it << " ";
                }
                logger.warn << std::endl;
                FATAL_ERROR_THROW;
        }
    }
}

template<>
const String &ParamReader::PatchParam::get(const String &name) const {
    if (dataStr_.find(name) == dataStr_.end()) {
        logger.warn << "ParamReader::PatchParam::get() - Cannot find \"" << name << "\"<String>" << std::endl;
        FATAL_ERROR_THROW;
    }
    return dataStr_.at(name);
}

template<>
const Scalar &ParamReader::PatchParam::get(const String &name) const {
    if (dataScl_.find(name) == dataScl_.end()) {
        logger.warn << "ParamReader::PatchParam::get() - Cannot find \"" << name << "\"<Scalar>" << std::endl;
        FATAL_ERROR_THROW;
    }
    return dataScl_.at(name);
}

template<>
const Vector &ParamReader::PatchParam::get(const String &name) const {
    if (dataVtr_.find(name) == dataVtr_.end()) {
        logger.warn << "ParamReader::PatchParam::get() - Cannot find \"" << name << "\"<Vector>" << std::endl;
        FATAL_ERROR_THROW;
    }
    return dataVtr_.at(name);
}
