#include "solver.h"

/// local function
PhysicalVar::MacroVars strvec_to_macro_vars(const string_vector &str_vec);

/// DUGKS_INCOMPRESSIBLE

using DUGKS_INCOMPRESSIBLE_CHECK_POINT = CheckPoint<DUGKS_INCOMPRESSIBLE>;

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::init_field(const PhysicalVar::MacroVars &_var) {
    for (auto &cell : solver.CELLS) {
        cell.init(_var);
    }
}

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    std::vector<PhysicalVar::MacroVars> data;
    for (int i = 0; i < reader.line_num(); i++) {
        auto each_line = split(reader[i]);
        if (each_line.empty()) continue;
        switch (read_case) {
            case 0:
                if (each_line[0] == "case=") {
                    if (each_line[1] != solver.config.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    continue;
                }
                if (each_line[0] == "mesh=") {
                    if (each_line[1] != solver.phy_mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    data.resize(solver.phy_mesh.cell_num());
                    continue;
                }
                if (each_line[0] == "step=") {
                    solver.step = stoi(each_line[1]);
                    std::stringstream ss;
                    ss << "Set step=" << solver.step << " for " << solver.config.solver;
                    highlight_println(ss.str());
                    continue;
                }
                if (each_line[0] == cell_mark) {
                    read_case = 1;
                    continue;
                }
                break;
            case 1:
                if (each_line[0] == ":end") break;
                solver.get_cell(stoi(each_line[0])).init(strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
}

TP_func void DUGKS_INCOMPRESSIBLE_CHECK_POINT::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
          "mesh= " << solver.phy_mesh.name << "\n" <<
          "step= " << solver.step << "\n\n" << cell_mark << "\n";
    for (auto &cell : solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
        << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
        << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z << "\t"
        << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z << "\n";
    }
    fp << end_mark;
}

/// CDUGKS_INCOMPRESSIBLE

using CDUGKS_INCOMPRESSIBLE_CHECK_POINT = CheckPoint<CDUGKS_INCOMPRESSIBLE>;

TP_func void CDUGKS_INCOMPRESSIBLE_CHECK_POINT::init_field(const PhysicalVar::MacroVars &_var) {
    for (auto &cell : solver.CELLS) {
        cell.init(_var);
    }
}

TP_func void CDUGKS_INCOMPRESSIBLE_CHECK_POINT::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    std::vector<PhysicalVar::MacroVars> data;
    for (int i = 0; i < reader.line_num(); i++) {
        auto each_line = split(reader[i]);
        if (each_line.empty()) continue;
        switch (read_case) {
            case 0:
                if (each_line[0] == "case=") {
                    if (each_line[1] != solver.config.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    continue;
                }
                if (each_line[0] == "mesh=") {
                    if (each_line[1] != solver.phy_mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    data.resize(solver.phy_mesh.cell_num());
                    continue;
                }
                if (each_line[0] == "step=") {
                    solver.step = stoi(each_line[1]);
                    std::stringstream ss;
                    ss << "Set step=" << solver.step << " for " << solver.config.solver;
                    highlight_println(ss.str());
                    continue;
                }
                if (each_line[0] == cell_mark) {
                    read_case = 1;
                    continue;
                }
                break;
            case 1:
                if (each_line[0] == ":end") break;
                solver.get_cell(stoi(each_line[0])).init(strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
}

TP_func void CDUGKS_INCOMPRESSIBLE_CHECK_POINT::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
       "mesh= " << solver.phy_mesh.name << "\n" <<
       "step= " << solver.step << "\n\n" << cell_mark << "\n";
    for (auto &cell : solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
           << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
           << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z << "\t"
           << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z << "\n";
    }
    fp << end_mark;
}

/// DUGKS_SHAKHOV

using DUGKS_SHAKHOV_CHECK_POINT = CheckPoint<DUGKS_SHAKHOV>;

TP_func void DUGKS_SHAKHOV_CHECK_POINT::init_field(const PhysicalVar::MacroVars &_var) {
    for (auto &cell : solver.CELLS) {
        cell.init(_var);
    }
}

TP_func void DUGKS_SHAKHOV_CHECK_POINT::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    std::vector<PhysicalVar::MacroVars> data;
    for (int i = 0; i < reader.line_num(); i++) {
        auto each_line = split(reader[i]);
        if (each_line.empty()) continue;
        switch (read_case) {
            case 0:
                if (each_line[0] == "case=") {
                    if (each_line[1] != solver.config.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    continue;
                }
                if (each_line[0] == "mesh=") {
                    if (each_line[1] != solver.phy_mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    data.resize(solver.phy_mesh.cell_num());
                    continue;
                }
                if (each_line[0] == "step=") {
                    solver.step = stoi(each_line[1]);
                    std::stringstream ss;
                    ss << "Set step=" << solver.step << " for " << solver.config.solver;
                    highlight_println(ss.str());
                    continue;
                }
                if (each_line[0] == cell_mark) {
                    read_case = 1;
                    continue;
                }
                break;
            case 1:
                if (each_line[0] == ":end") break;
                solver.get_cell(stoi(each_line[0])).init(strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
}

TP_func void DUGKS_SHAKHOV_CHECK_POINT::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
       "mesh= " << solver.phy_mesh.name << "\n" <<
       "step= " << solver.step << "\n\n" << cell_mark << "\n";
    for (auto &cell : solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
           << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
           << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z << "\t"
           << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z << "\n";
    }
    fp << end_mark;
}


/// local function
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
