#include "main.h"

int main(int argc, char** argv) {
    meso_init();

    ArgParser parser(argc, argv);
    debug_mode = parser.parse_switch("debug");

    std::string parsed_string;

    /// mesh
    parsed_string = parser.parse_param<std::string>("parse_mesh", STRING_NULL);
    if (parsed_string != STRING_NULL) return handle_parse_mesh(parsed_string);

    /// case
    parsed_string = parser.parse_param<std::string>("case", STRING_NULL);
    if (parsed_string != STRING_NULL) {
        ConfigReader config("./config/" + parsed_string);
        if (config.solver == "dugks@incompressible") return handle_solver<DUGKS_INCOMPRESSIBLE>(config, parser);
        if (config.solver == "dugks@shakhov") return handle_solver<DUGKS_SHAKHOV>(config, parser);
    }

    if (debug_mode) debug_println("exit");

    return 0;
}

int handle_parse_mesh(const std::string &path) {
    MESH::ListMesh mesh(MESH_TYPE_NORMAL, "parsed_mesh");
    mesh.load(path);
    mesh.build();
    mesh.info();
    /// output
    MeshWriter<MESH::ListMesh> writer{"./mesh.dat", mesh};
    writer.write_head({"volume"});
    writer.write_node();
    std::vector<double> data(mesh.cell_num());
    for (int i = 0; i < mesh.cell_num(); i++) data[i] = mesh.CELLS[i].volume;
    writer.write_data(data);
    writer.write_geom();
    writer.close();
    return 0;
}
