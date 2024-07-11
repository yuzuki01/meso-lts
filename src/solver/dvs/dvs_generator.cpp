#include "solver/solver.h"


using namespace MESO;

StringList read_dvs_lines(const std::string &file_path) {
    StringList lines;
    {
        MESO::FileReader::BasicReader reader(file_path);
        lines = reader.read_lines();
    }
    return lines;
}

template <>
fvmMesh::Mesh DVS::generate_dvs(const std::string &file_path, GaussHermiteParams &params) {
    auto lines = read_dvs_lines(file_path);
    const int n = static_cast<int>(lines.size()) - 1;
    auto data = MESO::Utils::split(lines[0]);
    if (data[0] != "Gauss-Hermite") {
        logger.error << "generate_dvs<Gauss-Hermite> get wrong type: " << data[0] << std::endl;
        throw std::invalid_argument("generate_dvs type error.");
    }
    /**
     * File format:
     *  1| Gauss-Hermite
     *  2| <x_0> <w_0>
     *  3| <x_1> <w_1>
     *  4| ...
     **/
    fvmMesh::Mesh mesh;
    const double x_scale = sqrt(2.0 * params.RT);
    const double w_scale = 1.0 / sqrt(M_PI);
    ScalarList x_list(n), w_list(n);
    for (int i = 0; i < n; ++i) {
        data = MESO::Utils::split(lines[i + 1]);    // skip first line
        x_list[i] = x_scale * std::stod(data[0]);
        w_list[i] = w_scale * std::stod(data[1]);
    }

    double max_mag = 0.0;
    if (params.dimension == 2) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                fvmMesh::Cell cell(i + j * n);
                cell.position = Vector{x_list[i], x_list[j], 0.0};
                cell.volume = w_list[i] * w_list[j];
                mesh.cells.push_back(cell);
                double cell_mag = cell.position.magnitude();
                if (max_mag < cell_mag) max_mag = cell_mag;
            }
        }
    } else {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    fvmMesh::Cell cell(i + (j + k * n) * n);
                    cell.position = Vector{x_list[i], x_list[j], x_list[k]};
                    cell.volume = x_list[i] * x_list[j] * x_list[k];
                    mesh.cells.push_back(cell);
                    double cell_mag = cell.position.magnitude();
                    if (max_mag < cell_mag) max_mag = cell_mag;
                }
            }
        }
    }
    mesh.update_num();
    mesh.max_cell_magnitude = max_mag;
    return mesh;
}

typedef std::unordered_map<int, std::vector<double>> ValueListMap;

ValueListMap nc_weight = {
        {1, {1. / 2.,    1. / 2.}},
        {2, {1. / 6.,    4. / 6.,   1. / 6.}},
        {3, {1. / 8.,    3. / 8.,   3. / 8.,    1. / 8.}},
        {4, {7. / 90.,   16. / 45., 2. / 15.,   16. / 45.,  7. / 90.}},
        {5, {19. / 288., 25. / 96., 25. / 144., 25. / 144., 25. / 96., 19. / 288.}},
        {6, {41. / 840., 9. / 35.,  9. / 280.,  34. / 105., 9. / 280., 9. / 35., 41. / 840.}}
};

template <>
fvmMesh::Mesh DVS::generate_dvs(const std::string &file_path, NewtonCotesParams &params) {
    auto lines = read_dvs_lines(file_path);
    auto data = MESO::Utils::split(lines[0]);
    if (data[0] != "Newton-Cotes") {
        logger.error << "generate_dvs<Newton-Cotes> get wrong type: " << data[0] << std::endl;
        throw std::invalid_argument("generate_dvs type error.");
    }
    /**
     * File format:
     *  1| Newton-Cotes
     *  2| <n> <mount> <scale>
     **/

    data = MESO::Utils::split(lines[1]);
    const int n = std::stoi(data[0]);
    const int mount = std::stoi(data[1]);
    const double scale = std::stod(data[2]);
    int total_point = n * mount + 1;
    auto &p_weight = nc_weight[n];
    double dh = 2.0 * scale / mount;
    ScalarList coordinate(total_point), weight(total_point);
    fvmMesh::Mesh mesh;
    for (int i = 0; i < total_point; i++) {
        coordinate[i] = scale * double(2 * i + 1 - total_point) / double(total_point - 1);
        weight[i] = 0.0;
    }
    for (int i = 0; i < mount; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            weight[i * n + j] += p_weight[j] * dh;
        }
    }
    double max_mag = 0.0;
    if (params.dimension == 2) {
        for (int j = 0; j < total_point; ++j) {
            for (int i = 0; i < total_point; ++i) {
                fvmMesh::Cell cell(i + j * total_point);
                cell.position = {coordinate[i], coordinate[j], 0.0};
                cell.volume = weight[i] * weight[j];
                mesh.cells.push_back(cell);
                double cell_mag = cell.position.magnitude();
                if (max_mag < cell_mag) max_mag = cell_mag;
            }
        }
    } else {
        for (int k = 0; k < total_point; ++k) {
            for (int j = 0; j < total_point; ++j) {
                for (int i = 0; i < total_point; ++i) {
                    fvmMesh::Cell cell(i + (j + k * total_point) * total_point);
                    cell.position = {coordinate[i], coordinate[j], coordinate[k]};
                    cell.volume = weight[i] * weight[j] * weight[k];
                    mesh.cells.push_back(cell);
                    double cell_mag = cell.position.magnitude();
                    if (max_mag < cell_mag) max_mag = cell_mag;
                }
            }
        }
    }
    mesh.update_num();
    mesh.max_cell_magnitude = max_mag;
    return mesh;
}
