#include "solver.h"

PhysicalVar::MacroVars strvec_to_macro_vars(const string_vector &str_vec) {
    PhysicalVar::MacroVars result;
    if (str_vec.size() < 9) {
        std::string err_str = "strvec_to_macro_vars length error";
        warn_println(err_str);
        throw std::invalid_argument(err_str);
    }
    result.density = stod(str_vec[1]);
    result.temperature = stod(str_vec[2]);
    result.velocity = {stod(str_vec[3]), stod(str_vec[4]), stod(str_vec[5])};
    result.heat_flux = {stod(str_vec[6]), stod(str_vec[7]), stod(str_vec[8])};
    return result;
}

/// DUGKS_INCOMPRESSIBLE

using DUGKS_INCOMPRESSIBLE_CHECK_POINT = CheckPoint<DUGKS_INCOMPRESSIBLE>;

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::init_field(const PhysicalVar::MacroVars &_var) {
    for (auto &cell : solver.CELLS) {
        cell.init(_var);
    }
    solver.logger << "init field from static/inlet var.";
    solver.logger.note();
}

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    for (int i = 0; i < reader.line_num(); i++) {
        string_vector each_line = split(reader[i]);
        if (each_line.empty() || std::strncmp(each_line[0].c_str(), "#", 1) == 0) continue;
        switch (read_case) {
            case 0:
                if (each_line[0] == "case=") {
                    if (each_line[1] != solver.config.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "mesh=") {
                    if (each_line[1] != solver.phy_mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "solver=") {
                    if (each_line[1] != solver.config.solver) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "step=") {
                    solver.step = stoi(each_line[1]);
                    std::stringstream ss;
                    ss << "  set step=" << solver.step << " for " << solver.config.solver;
                    info_println(ss.str());
                } else if (each_line[0] == "cell:") {
                    read_case = 1;
                }
                break;
            case 1:
                if (each_line[0] == ":end") {
                    read_case = 0;
                    break;
                }
                solver.get_cell(stod(each_line[0])).init(strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
    solver.logger << "init field from file: " << file_path;
    solver.logger.note();
}

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
       "mesh= " << solver.phy_mesh.name << "\n" <<
       "solver= " << solver.config.solver << "\n" <<
       "step= " << solver.step << "\n\ncell:\n";
    for (auto &cell : solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
           << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
           << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z
           << "\t"
           << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z
           << "\n";
    }
    fp << ":end";
}


/// DUGKS_SHAKHOV

using DUGKS_SHAKHOV_CHECK_POINT = CheckPoint<DUGKS_SHAKHOV>;

TP_func void DUGKS_SHAKHOV_CHECK_POINT::init_field(const PhysicalVar::MacroVars &_var) {
    for (auto &cell : solver.CELLS) {
        cell.init(_var);
    }
    solver.logger << "init field from static/inlet var.";
    solver.logger.note();
}

TP_func void DUGKS_SHAKHOV_CHECK_POINT::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    for (int i = 0; i < reader.line_num(); i++) {
        string_vector each_line = split(reader[i]);
        if (each_line.empty() || std::strncmp(each_line[0].c_str(), "#", 1) == 0) continue;
        switch (read_case) {
            case 0:
                if (each_line[0] == "case=") {
                    if (each_line[1] != solver.config.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "mesh=") {
                    if (each_line[1] != solver.phy_mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "solver=") {
                    if (each_line[1] != solver.config.solver) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                } else if (each_line[0] == "step=") {
                    solver.step = stoi(each_line[1]);
                    std::stringstream ss;
                    ss << "  set step=" << solver.step << " for " << solver.config.solver;
                    info_println(ss.str());
                } else if (each_line[0] == "cell:") {
                    read_case = 1;
                }
                break;
            case 1:
                if (each_line[0] == ":end") {
                    read_case = 0;
                    break;
                }
                solver.get_cell(stod(each_line[0])).init(strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
    solver.logger << "init field from file: " << file_path;
    solver.logger.note();
}

TP_func void DUGKS_SHAKHOV_CHECK_POINT::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
       "mesh= " << solver.phy_mesh.name << "\n" <<
       "solver= " << solver.config.solver << "\n" <<
       "step= " << solver.step << "\n\ncell:\n";
    for (auto &cell : solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
           << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
           << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z
           << "\t"
           << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z
           << "\n";
    }
    fp << ":end";
}
