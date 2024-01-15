#include "solver.h"

#ifdef SOLVER_GKS

using Scheme = GKS;
using SCell = Scheme::Cell;
using SFace = Scheme::Face;


GKS::GKS(ConfigReader &_config, ArgParser &_parser) : BasicSolver(_config, _parser), check_point(*this) {
    /// Read config
    Rho0 = config.get<double>("density");
    Ma = config.get<double>("Ma");
    R = config.get<double>("R");
    T0 = config.get<double>("temperature");
    L = config.get<double>("length");
    CFL = config.get<double>("CFL");
    K = config.get<int>("K");
    /// Mesh
    mesh.load("./mesh/" + config.phy_mesh);
    mesh.build();
    D = mesh.dimension();
}

void Scheme::info() const {
    mesh.info();
}

GKS::Cell::Cell(MESH::Cell &cell, GKS::Scheme &_solver) : mesh_cell(cell), solver(_solver) {}

GKS::Face::Face(MESH::Face &face, GKS::Scheme &_solver) : mesh_face(face), solver(_solver) {}


/// Physical Formula
double Scheme::lambda(double temperature) const {
    /// lambda = rho / (2 * p) = 1 / (2 * R * T)
    return 1.0 / (2.0 * R * temperature);
}

void SCell::init(const Physical::MacroVars &init_var) {
    lsp = generate_least_square(mesh_cell.key, solver.mesh);
}

void SCell::update() {
    ///todo
}

void SFace::reconstruct() {
    if (mesh_face.boundary_type != MeshBC_interface) return boundary();
    /**
     * The 1st Order KFVS (Kinetic Flux Vector Splitting)
     **/
    auto &left_cell = solver.get_cell(mesh_face.on_cell_key);
    auto &right_cell = solver.get_cell(mesh_face.inv_cell_key);
    auto &nv = mesh_face.on_cell_nv;
    /// cell and vector follows as:
    /// | left_cell |*| right_cell |
    ///   |*| is where the face is now
    /// the positive direction is pointing to the right cell (->) and along nv in real mesh

    /// coordinate transfer (global -> local)
    ///     M * u = {|u|, 0, 0}
    Vec3D U_l, U_r;
    auto &M = mesh_face.rotate_matrix;
    auto &Mi = mesh_face.inv_rotate_matrix;
    U_l = M * left_cell.macro_vars.velocity;
    U_r = -(M * right_cell.macro_vars.velocity);
    /// KFVS
    double rho_l = left_cell.macro_vars.density, rho_r = left_cell.macro_vars.density;
    double lambda_l = solver.lambda(left_cell.macro_vars.temperature),
            lambda_r = solver.lambda(right_cell.macro_vars.temperature);
    /// todo
}

void SFace::boundary() {
    switch (mesh_face.boundary_type) {
        case MeshBC_isothermal_wall: {

        }
            break;
        default: {
            solver.logger << "Caught unsupported boundary type: " << mesh_face.info();
            solver.logger.error();
            throw std::invalid_argument("boundary error");
        }
    }
}

void Scheme::do_step() {
#pragma omp parallel for
    for (int i = 0; i < int(mesh.face_num()); i++) {
        auto &face = get_face(i);
        face.reconstruct();
    }

#pragma omp parallel for
    for (int i = 0; i < int(mesh.cell_num()); ++i) {
        auto &cell = get_cell(i);
        cell.update();
    }

    /// OpenMP loop end
    step++;
    if (step >= max_step && !is_crashed) {
        continue_to_run = false;
        logger << "reach max_step=" << max_step;
        logger.highlight();
        do_residual();
        do_save();
        return;
    }
    if (step % save_interval == 0) {
        do_residual();
        do_save();
    }
}

void Scheme::init() {
    /// if no init, then cannot run
    continue_to_run = true;
    /// Compute Params
    gamma = (5.0 + K) / (3.0 + K);
    double sound = sqrt(gamma * R * T0);
    dt = CFL * mesh.min_mesh_size / ((Ma > 1.0 ? Ma : 1.0) * sound);
    /// Boundary
    if (!config.set_mesh_mark(mesh)) {
        continue_to_run = false;
        return;
    }
    logger << "solver info:";
    logger.note();
    data_double_println({"density", "temperature", "length", "R"},
                        {Rho0, T0, L, R});
    data_double_println({"CFL", "min-cell-size", "time-step"},
                        {CFL, mesh.min_mesh_size, dt});

    /// Create Scheme Objects
    logger << "Create scheme-objects.";
    logger.note();
    for (auto &it: mesh.CELLS) {
        CELLS.emplace_back(it, *this);
    }

    auto check_point_file = parser.parse_param<std::string>("check_point", STRING_NULL);
    if (check_point_file == STRING_NULL) {
        Physical::MacroVars init_var{};
        init_var.density = Rho0;
        init_var.temperature = T0;
        init_var.velocity = {0.0, 0.0, 0.0};
        for (auto &mark: mesh.MARKS) {
            if (MESH::MarkTypeID[mark.type] == MeshBC_inlet) {
                init_var.density = mark.density;
                init_var.temperature = mark.temperature;
                init_var.velocity = mark.velocity;
                break;
            }
        }
        for (auto &mark: mesh.MARKS) {
            if (MESH::MarkTypeID[mark.type] == MeshBC_inlet) {
                init_var.density = mark.density;
                init_var.temperature = mark.temperature;
                init_var.velocity = mark.velocity;
                break;
            }
        }
        check_point.init_field(init_var);
    } else {
        check_point.init_from_file(check_point_file);
    }

    logger << "    scheme-objects: cell - ok.";
    logger.info();
    for (auto &face: mesh.FACES) {
        FACES.emplace_back(face, *this);
        /// face init
    }
    logger << "    scheme-objects: face - ok.";
    logger.info();

    omp_set_num_threads(config.thread);
    logger << "init - ok.";
    logger.note();
}


void Scheme::do_save() {
    std::stringstream ss;
    ss << "./result/" << config.name << "/step_" << step << ".dat";
    string_vector values;
    if (mesh.dimension() == 2) values = {"density", "temperature", "velocity-x", "velocity-y"};
    else values = {"density", "temperature", "velocity-x", "velocity-y", "velocity-z"};
    MeshWriter<MESH::StaticMesh> writer(ss.str(), mesh);
    writer.write_head(values);
    writer.write_node();

    std::vector<double> data(mesh.cell_num());

    /// density
    for (auto &it: CELLS) data[it.mesh_cell.key] = it.macro_vars.density;
    writer.write_data(data);
    data.clear();
    /// temperature
    for (auto &it: CELLS) data[it.mesh_cell.key] = it.macro_vars.temperature;
    writer.write_data(data);
    data.clear();
    /// velocity-x
    for (auto &it: CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.x;
    writer.write_data(data);
    data.clear();
    /// velocity-y
    for (auto &it: CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.y;
    writer.write_data(data);
    data.clear();
    /// velocity-z
    if (D == 3) {
        for (auto &it: CELLS) data[it.mesh_cell.key] = it.macro_vars.velocity.z;
        writer.write_data(data);
        data.clear();
    }

    writer.write_geom();
    writer.close();

    /// check point
    std::stringstream ss2;
    ss2 << "./result/" << config.name << "/" << config.name << ".check_point";
    check_point.write_to_file(ss2.str());

    logger << "Save to file: " << ss.str();
    logger.highlight();
}

SCell &Scheme::get_cell(const int &_key) {
    return CELLS[_key];
}

SFace &Scheme::get_face(const int &_key) {
    return FACES[_key];
}


void Scheme::do_residual() {
    std::vector<double> residual(2);
    double sum1, sum2;
    /// density
    sum1 = sum2 = 0.0;
    for (auto &cell: CELLS) {
        double a, b;
        a = cell.macro_vars.density * cell.mesh_cell.volume;
        b = cell.macro_vars.old_density * cell.mesh_cell.volume;
        sum1 += fabs((a - b) * (a - b));
        sum2 += fabs(b * b);
        cell.macro_vars.old_density = cell.macro_vars.density;
    }
    residual[0] = sqrt(sum1 / sum2);
    /// energy
    sum1 = sum2 = 0.0;
    for (auto &cell: CELLS) {
        double a, b;
        a = cell.macro_vars.energy * cell.mesh_cell.volume;
        b = cell.macro_vars.old_energy * cell.mesh_cell.volume;
        sum1 += fabs((a - b) * (a - b));
        sum2 += fabs(b * b);
        cell.macro_vars.old_energy = cell.macro_vars.energy;
    }
    residual[1] = sqrt(sum1 / sum2);
    logger << "Residual at step=" << step;
    logger.info();
    data_sci_double_println({"density", "energy"}, residual);
}

/// check_point

TP_func void CheckPoint<Scheme>::init_field(const Physical::MacroVars &_var) {
    for (auto &cell: solver.CELLS) {
        cell.init(_var);
    }
}

TP_func void CheckPoint<Scheme>::init_from_file(const std::string &file_path) {
    BasicReader reader("check_point", file_path);
    int read_case = 0;
    std::vector<Physical::MacroVars> data;
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
                    if (each_line[1] != solver.mesh.name) {
                        solver.continue_to_run = false;
                        warn_println("check_point error.");
                        return;
                    }
                    data.resize(solver.mesh.cell_num());
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
                solver.get_cell(stoi(each_line[0])).init(Physical::strvec_to_macro_vars(each_line));
                break;
            default:
                break;
        }
    }
}

TP_func void CheckPoint<Scheme>::write_to_file(const std::string &file_path) {
    std::ofstream fp;
    fp.open(file_path, std::ios::out | std::ios::trunc);
    fp << "case= " << solver.config.name << "\n" <<
       "mesh= " << solver.mesh.name << "\n" <<
       "step= " << solver.step << "\n\n" << cell_mark << "\n";
    for (auto &cell: solver.CELLS) {
        fp << cell.mesh_cell.key << "\t" << std::setprecision(DATA_PRECISION)
           << cell.macro_vars.density << "\t" << cell.macro_vars.temperature << "\t"
           << cell.macro_vars.velocity.x << "\t" << cell.macro_vars.velocity.y << "\t" << cell.macro_vars.velocity.z
           << "\t"
           << cell.macro_vars.heat_flux.x << "\t" << cell.macro_vars.heat_flux.y << "\t" << cell.macro_vars.heat_flux.z
           << "\n";
    }
    fp << end_mark;
}

/// MMDF

GKS::MMDF::MMDF(double _u, double _lambda, int _k, int integral_zone) : u(_u), lambda(_lambda), K(_k),
zone(integral_zone) {
    for (int i = 0; i < 7; ++i) {
        is_computed[i] = false;
        _moment[i] = 0.0;
    }
}

double GKS::MMDF::operator()(int _o) {
    switch (_o) {
        case 0:
            if (!is_computed[_o]) {
                switch (zone) {
                    case full_zone:
                        _moment[_o] = 1.0;
                        break;
                    case left_zone:
                        _moment[_o] = erfc(-sqrt(lambda * u)) / 2.0;
                        break;
                    case right_zone:
                        _moment[_o] = erfc(sqrt(lambda) * u) / 2.0;
                        break;
                }
                is_computed[_o] = true;
            }
            break;
        case 1:
            if (!is_computed[_o]) {
                switch (zone) {
                    case full_zone:
                        _moment[_o] = u;
                        break;
                    case left_zone:
                        _moment[_o] = u * operator()(0) + exp(-lambda * u * u) / (2.0 * sqrt(Pi * lambda));
                        break;
                    case right_zone:
                        _moment[_o] = u * operator()(0) - exp(-lambda * u * u) / (2.0 * sqrt(Pi * lambda));
                        break;
                }
                is_computed[_o] = true;
            }
            break;
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
            if (!is_computed[_o]) {
                _moment[_o] = u * operator()(_o - 1) + (_o - 1.0) / (2.0 * lambda) * operator()(_o - 2);
                is_computed[_o] = true;
            }
            break;
        default:
            warn_println("MMDF error.");
            return 0.0;
    }
    return _moment[_o];
}

double GKS::MMDF::operator[](int _o) const {
    switch (_o) {
        case 0:
            return 1.0;
        case 2:
            return K / (2.0 * lambda);
        default:
            return 0.0;
    }
}

#endif