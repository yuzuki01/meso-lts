#include "solver.h"

ConfigReader::ConfigReader(const std::string &file_path) : BasicReader("ConfigReader", file_path) {
    // parse context
    int read_case = 0;
    for (int i = 0; i < int(lines.size()); i++) {
        auto data = split(lines[i]);
        /// skip empty
        if (data.empty() || std::strncmp(data[0].c_str(), "#", 1) == 0) continue;
        switch (read_case) {
            case 0:
                if (data[0] == "case-info:") read_case = 1;
                if (data[0] == def_mark) read_case = 2;
                if (data[0] == bc_mark) read_case = 3;
                break;
            case 1: {
                int loop = 0;
                while (true) {
                    data = split(lines[i]);
                    if (data.empty() || std::strncmp(data[0].c_str(), "#", 1) == 0) {
                        loop++;
                        if (loop >= 1000) {
                            warn_println("Reader caught too many empty lines.");
                            break;  // break from while
                        }
                        i++;
                        continue;
                    }
                    if (data[0] == end_mark) break;  // break from while
                    else if (data[0] == "name") name = data[1];
                    else if (data[0] == "phy_mesh") phy_mesh = data[1];
                    else if (data[0] == "dvs_mesh") dvs_mesh = data[1];
                    else if (data[0] == "solver") solver = data[1];
                    else if (data[0] == "thread") thread = stoi(data[1]);
                    i++;
                }
                read_case = 0;
                break;  // break from switch-config
            }
            case 2: {
                int loop = 0;
                while (true) {
                    data = split(lines[i]);
                    if (data.empty() || std::strncmp(data[0].c_str(), "#", 1) == 0) {
                        loop++;
                        if (loop >= 1000) {
                            warn_println("Reader caught too many empty lines.");
                            break;  // break from while
                        }
                        i++;
                        continue;
                    }
                    if (data[0] == end_mark) break;  // break from while
                    param[data[0]] = data[1];
                    i++;
                }
                read_case = 0;
                break;  // break from switch-config
            }
            case 3: {
                int loop = 0;
                BoundaryParam bc_param;
                while (true) {
                    data = split(lines[i]);
                    if (data.empty() || std::strncmp(data[0].c_str(), "#", 1) == 0) {
                        loop++;
                        if (loop >= 1000) {
                            warn_println("Reader caught too many empty lines.");
                            break;  // break from while
                        }
                        i++;
                        continue;
                    }
                    if (data[0] == end_mark) break;  // break from while
                    else if (data[0] == "name") bc_param.name = data[1];
                    else if (data[0] == "type") bc_param.type = data[1];
                    else if (data[0] == "density") bc_param.density = stod(data[1]);
                    else if (data[0] == "temperature") bc_param.temperature = stod(data[1]);
                    else if (data[0] == "velocity-x") bc_param.velocity.x = stod(data[1]);
                    else if (data[0] == "velocity-y") bc_param.velocity.y = stod(data[1]);
                    else if (data[0] == "velocity-z") bc_param.velocity.z = stod(data[1]);
                    i++;
                }
                mark_sets.push_back(bc_param);
                std::stringstream ss;
                ss << "Read mark<" << bc_param.name << "> params from config file.";
                info_println(ss.str());
                read_case = 0;
                break;  // break from switch-config
            }
            default:
                break;
        }
    }
    logger("Case-" + name);
    mark_sets.shrink_to_fit();
    info_println("Read config: " + file_path + " - ok.");
}

void ConfigReader::info() {
    logger << "Config info:";
    logger.note();
    std::stringstream ss;
    ss << "  config name: " << name << "\n";
    ss << "   phy mesh: " << phy_mesh << "\n";
    ss << "   dvs mesh: " << dvs_mesh << "\n";
    ss << "    solver : " << solver << "\n";
    ss << "    thread : " << thread << "\n";
    info_println(ss.str());
}

std::string ConfigReader::operator[](const std::string &key) {
    auto it = param.find(key);
    if (it != param.end()) return param[key];
    return STRING_NULL;
}

template<> std::string ConfigReader::get<std::string>(const std::string &_key) {
    return this->operator[](_key);
}

template<> double ConfigReader::get<double>(const std::string &_key) {
    std::string s = this->operator[](_key);
    if (s == STRING_NULL) {
        warn_println("ConfigReader cannot find key: " + _key + "<double>, return 0.0");
        return 0.0;
    }
    return stod(s);
}

template<> int ConfigReader::get<int>(const std::string &_key) {
    std::string s = this->operator[](_key);
    if (s == STRING_NULL) {
        warn_println("ConfigReader cannot find key: " + _key + "<int>, return 0");
        return 0;
    }
    return stoi(s);
}

template<> bool ConfigReader::get<bool>(const std::string &_key) {
    std::string s = this->operator[](_key);
    if (s == STRING_NULL) {
        warn_println("ConfigReader cannot find key: " + _key + "<int>, return 0");
        return false;
    }
    if (s == "true" || s == "True") return true;
    return false;
}

bool ConfigReader::set_mesh_mark(MESH::StaticMesh &mesh) {
    bool flag = true;
    for (auto &bc : mark_sets) {
        if (!mesh.set_mark(bc)) flag = false;
    }
    return flag;
}

void ConfigReader::set_mesh_mark(MESH::MapMesh &mesh) {
    for (auto &bc : mark_sets) mesh.set_mark(bc);
}
