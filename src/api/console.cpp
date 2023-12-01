#include "api.h"

int handle_parse_mesh(const std::string &path) {
    MESH::StaticMesh mesh(MESH_TYPE_NORMAL, "parsed_mesh");
    mesh.load(path);
    mesh.build();
    mesh.info();
    /// output
    MeshWriter<MESH::StaticMesh> writer{"./mesh.dat", mesh};
    writer.write_head({"volume"});
    writer.write_node();
    std::vector<double> data(mesh.cell_num());
    for (int i = 0; i < mesh.cell_num(); i++) data[i] = mesh.CELLS[i].volume;
    writer.write_data(data);
    writer.write_geom();
    writer.close();
    return 0;
}

int handle_help() {
    Logger logger("meso");
    logger << " meso -<switch> --<param> <value>" << "\n" <<
    "     switch:" << "\n" <<
    "       h / help        -- help information" << "\n" <<
    "       debug           -- debug symbol" << "\n" <<
    "     param:" << "\n" <<
    "       case            -- fetch config file in ./config/" << "\n" <<
    "       max_step        -- set max step" << "\n" <<
    "       save_interval   -- set save step interval" << "\n" <<
    "       parse_mesh      -- output tecplot file of mesh file" << "\n" <<
    "       ...";
    logger.info();
    return 0;
}
