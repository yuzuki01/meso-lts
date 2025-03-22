#include "solver/solver.h"

using namespace MESO::Solver;

const std::unordered_map<std::string, int> MESO::Solver::mark_type_map = {
        {"fluid-interior",   BoundaryType::fluid_interior},
        {"farfield-inlet",   BoundaryType::farfield_inlet},
        {"pressure-inlet",   BoundaryType::pressure_inlet},
        {"freestream-inlet", BoundaryType::freestream_inlet},
        {"farfield-outlet",  BoundaryType::farfield_outlet},
        {"pressure-outlet",  BoundaryType::pressure_outlet},
        {"wall",             BoundaryType::wall},
        {"slip-wall",        BoundaryType::slip_wall},
        {"symmetry",         BoundaryType::symmetry},
};


Config::Config(const std::string &file_path) : file_path(file_path) {
    {
        auto mark_ptr = new Mark;
        mark_ptr->id = 1;
        mark_ptr->name = mark_ptr->type_name = "fluid-interior";
        mark_ptr->type = mark_type_map.at(mark_ptr->type_name);
        marks[mark_ptr->name] = *mark_ptr;
    }
    update_config(true);
    logger.info << "Loaded case file: " << file_path << std::endl;
}

void Config::update_config(bool update_patch) {
    MESO::FileReader::BasicReader reader(file_path);
    StringList lines = reader.read_lines();
    StringList data;
    const int line_size = int(lines.size());
    for (int i = 0; i < line_size; ++i) {
        data = MESO::Utils::split(lines[i]);
        if (data[0] == "#") continue;       // skip #
        if (data[0] == "[settings]") {
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "#") {   // skip #
                    ++i;
                    continue;
                }
                if (data[0] == "[mark]" or data[0] == "[group]") {
                    --i;
                    break;
                }
                settings[data[0]] = data[1];
                ++i;
            }
            continue;
        } else if (data[0] == "[mark]") {
            Mark mark;
            mark.id = int(marks.size()) + 1;
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "#") {   // skip #
                    ++i;
                    continue;
                }
                if (data[0] == "[mark]" or data[0] == "[group]") {
                    --i;
                    break;
                }
                if (not update_patch) {
                    ++i;
                    continue;
                }
                if (data[0] == "name") mark.name = data[1];
                else if (data[0] == "type") {
                    mark.type_name = data[1];
                    mark.type = mark_type_map.at(mark.type_name);
                } else mark.patch.set_value(data);
                ++i;
            }
            marks[mark.name] = mark;
            continue;
        } else if (data[0] == "[group]") {
            Group group;
            group.id = int(groups.size()) + 1;
            ++i;
            while (i < line_size) {
                data = MESO::Utils::split(lines[i]);
                if (data[0] == "#") {   // skip #
                    ++i;
                    continue;
                }
                if (data[0] == "[mark]" or data[0] == "[group]") {
                    --i;
                    break;
                }
                if (not update_patch) {
                    ++i;
                    continue;
                }
                if (data[0] == "name") {
                    group.name = data[1];
                } else group.patch.set_value(data);
                ++i;
            }
            groups[group.name] = group;
            continue;
        }
    }
}

template<>
std::string Config::get(const ConfigKey &key, std::string default_value, bool print_error) {
    if (settings.find(key) == settings.end()) {
        if (print_error)
            logger.warn << "Cannot found settings[" << key << "], use default value instead. " << default_value
                        << std::endl;
        return default_value;
    }
    return settings[key];
}

template<>
int Config::get(const ConfigKey &key, int default_value, bool print_error) {
    if (settings.find(key) == settings.end()) {
        if (print_error)
            logger.warn << "Cannot found settings[" << key << "], use default value instead. " << default_value
                        << std::endl;
        return default_value;
    }
    return stoi(settings[key]);
}

template<>
double Config::get(const ConfigKey &key, double default_value, bool print_error) {
    if (settings.find(key) == settings.end()) {
        if (print_error)
            logger.warn << "Cannot found settings[" << key << "], use default value instead. " << default_value
                        << std::endl;
        return default_value;
    }
    return stod(settings[key]);
}

template<>
bool Config::get(const ConfigKey &key, bool default_value, bool print_error) {
    if (settings.find(key) == settings.end()) {
        if (print_error)
            logger.warn << "Cannot found settings[" << key << "], use default value instead. " << default_value
                        << std::endl;
        return default_value;
    }
    return (settings[key] == "True");
}

Group &Config::get_cell_group(MESO::fvmMesh::Cell &cell, MESO::fvmMesh::Mesh &mesh) {
    auto &key = mesh.cell_names[cell.group_id];
    if (groups.find(key) == groups.end()) {
        std::stringstream error_massage;
        error_massage << "get_cell_group() missing key=" << key;
        logger.warn << "[Error] " << error_massage.str() << std::endl;
        throw std::invalid_argument(error_massage.str());
    }
    return groups[key];
}

Mark &Config::get_face_group(MESO::fvmMesh::Face &face, MESO::fvmMesh::Mesh &mesh) {
    auto &key = mesh.face_names[face.group_id];
    if (marks.find(key) == marks.end()) {
        std::stringstream error_massage;
        error_massage << "get_face_group() missing key=" << key;
#pragma omp critical
        {
            logger.warn << "[Error] " << error_massage.str() << std::endl;
        }
        throw std::invalid_argument(error_massage.str());
    }
    return marks[key];
}

void Config::info() {
    logger.note << "Settings:\n";
    for (auto &[name, value]: settings) {
        logger.info << "  " << name << std::setw(10) << "\t\t" << value << "\n";
    }
    logger.note << "\nCell groups:\n";
    for (auto &[name, group]: groups) {
        logger.info << "  <" << name << ">\tid: " << group.id << "\n";
    }
    logger.note << "\nFace groups:\n";
    for (auto &[name, mark]: marks) {
        logger.info << "  <" << name << ">\n"
                    << "    id: " << mark.id << "\ttype: " << mark.type_name << " <" << mark.type << ">" << "\n";
    }
    logger.info << std::endl;
}
