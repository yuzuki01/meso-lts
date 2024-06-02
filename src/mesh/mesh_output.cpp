#include "mesh/mesh.h"

/// local func
void write_mesh_geom(std::fstream &fp, MESO::Mesh::Zone *mesh_ptr);

/// global func
using namespace MESO::Mesh;


void Zone::output(const std::string &file_name) {
    const int dim = dimension();
    const int DATA_PRECISION = 15;
    const int LINE_DATA_NUM = 30;
    std::fstream fp;
    std::stringstream ss;
    ss << file_name << ".msh.plt";
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << ss.str() << std::endl;
        fp.close();
        return;
    }
    // write head
    fp << R"(VARIABLES = "X", "Y")";
    if (dim == 3) {
        fp << R"( ,"Z")";
    }
    fp << R"(, "GroupId")";
    fp << std::endl;
    fp << "ZONE T=\"" << "MeshData" << "\", N=" << NNODE << ", E=" << NCELL << ", VARLOCATION=([1-" << dim
       << "]=NODAL, [3]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE="
       << (dim == 2 ? "FEQUADRILATERAL" : "FEBRICK")
       << std::endl;
    int count;
    // write node
    count = 0;
    fp << std::endl << "## Node-x" << std::endl;
    for (auto &node: nodes) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    count = 0;
    fp << std::endl << "## Node-y" << std::endl;
    for (auto &node: nodes) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    if (dim == 3) {
        count = 0;
        fp << std::endl << "## Node-z" << std::endl;
        for (auto &node: nodes) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
        fp << std::endl;
    }

    count = 0;
    fp << std::endl << "## GroupId" << std::endl;
    for (auto &cell: cells) {
        fp << "\t" << cell.group_id;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    // write geom
    write_mesh_geom(fp, this);

    // close
    fp.close();
    logger.note << "Write mesh to file: ";
    logger.info << ss.str() << std::endl;
}

void Zone::output(const std::string &file_name,
                  std::initializer_list<std::string> names,
                  std::initializer_list<MESO::Field<Scalar> *> values,
                  int step, double solution_time) {
    if (names.size() != values.size())
        throw std::invalid_argument("Wrong size when output result");
    const int vc = int(names.size());
    const int dim = dimension();
    const int DATA_PRECISION = 15;
    const int LINE_DATA_NUM = 30;
    std::fstream fp;
    std::stringstream ss;
    ss << file_name << ".dat.plt";
    fp.open(ss.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << ss.str() << std::endl;
        fp.close();
        return;
    }
    // write head
    fp << R"(VARIABLES = "X", "Y")";
    if (dim == 3) {
        fp << R"( ,"Z")";
    }
    for (int i = 0; i < vc; ++i) {
        fp << " ,\"" << *(names.begin() + i) << "\"";
    }
    fp << std::endl;

    fp << "ZONE T=\""
       << ((step >= 0) ? "SolutionData(step=" + std::to_string(step) + ")" : "SolutionData")
       << "\", N=" << NNODE << ", E=" << NCELL << ", VARLOCATION=([1-" << dim << "]=NODAL, [";
    if (vc == 1) {
        fp << dim + 1;
    } else {
        fp << dim + 1 << "-" << dim + vc;
    }
    fp << "]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE="
       << (dim == 2 ? "FEQUADRILATERAL" : "FEBRICK");
    if (solution_time >= 0.0) fp << ", SOLUTIONTIME=" << solution_time;
    fp << std::endl;
    int count;
    // write node
    count = 0;
    fp << std::endl << "## Node-x" << std::endl;
    for (auto &node: nodes) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    count = 0;
    fp << std::endl << "## Node-y" << std::endl;
    for (auto &node: nodes) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    if (dim == 3) {
        count = 0;
        fp << std::endl << "## Node-z" << std::endl;
        for (auto &node: nodes) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.position.z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
        fp << std::endl;
    }

    for (int i = 0; i < vc; ++i) {
        auto data_ptr = *(values.begin() + i);
        auto &data = *data_ptr;
        count = 0;
        fp << std::endl << "## CellVar-" << *(names.begin() + i) << std::endl;
        for (auto value: data.values) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << value;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
        fp << std::endl;
    }

    // write geom
    write_mesh_geom(fp, this);

    // close
    fp.close();
    logger.note << "Save result to file: ";
    logger.info << ss.str() << std::endl;
}

void write_mesh_geom(std::fstream &fp, MESO::Mesh::Zone *mesh_ptr) {
    // write geom
    fp << std::endl << "## geom" << std::endl;
    for (auto &cell: mesh_ptr->cells) {
        const int node_num = MESO::Geom::node_num(cell.geom_type);
        MESO::ObjectIdList node(node_num);
        for (int i = 0; i < node_num; ++i) {
            node[i] = cell.node_id[i] + 1;
        }
        switch (cell.geom_type) {
            case MESO::Geom::Quad:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[3];
                break;
            case MESO::Geom::Tria:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[0];
                break;
            case MESO::Geom::Brick:
                fp << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
                   << " " << node[4] << " " << node[5] << " " << node[7] << " " << node[6];
                break;
            case MESO::Geom::Wedge:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
                   << " " << node[3] << " " << node[4] << " " << node[5] << " " << node[5];
                break;
            case MESO::Geom::Tetra:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
                   << " " << node[3] << " " << node[3] << " " << node[3] << " " << node[3];
                break;
            case MESO::Geom::Pyram:
                fp << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
                   << " " << node[4] << " " << node[4] << " " << node[4] << " " << node[4];
                break;
            default:
                fp << " unsupported type.";
        }
        fp << std::endl;
    }
}
