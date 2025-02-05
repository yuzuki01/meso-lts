#include "mesh/mesh.h"


namespaceMesoMesh


void writeHead(std::fstream &fp, const fvMesh& mesh, const List<String>& varName);
void writeNode(std::fstream &fp, const fvMesh& mesh);
void writeGeom(std::fstream &fp, const fvMesh& mesh);

/**
 * ===================================================
 * --------------------- fvMesh ----------------------
 * ===================================================
 **/

void fvMesh::output() {
    if (MPI::rank != MPI::mainRank) return;
    Utils::mkdir("constant");
    const int dim = dimension();
    const int DATA_PRECISION = 15;
    const int LINE_DATA_NUM = 30;
    std::fstream fp;
    std::stringstream ss;
    fp.open("constant/grid.plt", std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << ss.str() << std::endl;
        fp.close();
        return;
    }
    // write head
    writeHead(fp, *this, {"Group", "Partition"});
    // write node
    writeNode(fp, *this);

    int count;
    List<ObjectId> cellValue(NCELL);

    for (int patchi = 0; patchi < zones_.size() ; ++patchi) {
        const auto& patch = zones_[patchi];
        for (const auto &ci : patch.group()) {
            cellValue[ci] = patchi;
        }
    }
    count = 0;
    fp << std::endl << "## Group" << std::endl;
    for (auto &cv: cellValue) {
        fp << "\t" << cv;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    for (int patchi = 0; patchi < parts_.size() ; ++patchi) {
        const auto& patch = parts_[patchi];
        for (const auto &ci : patch.group()) {
            cellValue[ci] = patchi;
        }
    }
    count = 0;
    fp << std::endl << "## Partition" << std::endl;
    for (auto &cv: cellValue) {
        fp << "\t" << cv;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    // write geom
    writeGeom(fp, *this);

    // close
    fp.close();

    /// Write cell/face/node files
    Utils::mkdir("constant/meshGeom");
    {
        List <Vector> node_pos;
        for (auto &it: nodes_) {
            node_pos.push_back(it.C());
        }
        Utils::output_list("constant/meshGeom/nodePosition", node_pos);
    }
    {
        List <Vector> cell_pos;
        List <Scalar> cell_volume;
        List <List<ObjectId>> cell_neighbor;
        for (auto &it: cells_) {
            cell_pos.push_back(it.C());
            cell_volume.push_back(it.V());
        }
        Utils::output_list("constant/meshGeom/cellPosition", cell_pos);
        Utils::output_list("constant/meshGeom/cellVolume", cell_volume);
    }
    {
        List <Vector> face_pos;
        List <Scalar> face_area;
        List <Set<ObjectId>> face_neighbor;
        for (auto &it: faces_) {
            face_pos.push_back(it.C());
            face_area.push_back(it.S());
        }
        Utils::output_list("constant/meshGeom/facePosition", face_pos);
        Utils::output_list("constant/meshGeom/faceArea", face_area);
    }

    logger.note << "Write mesh to file: ";
    logger.info << ss.str() << std::endl;
}


/**
 * ===================================================
 * ---------------------- local ----------------------
 * ===================================================
 **/

void writeHead(std::fstream &fp, const fvMesh& mesh, const List<String>& varName) {
    const Label dim = mesh.dimension();
    fp << "TITLE = \"Grid\"\n"
          "FileType = \"GRID\"\n";
    if (dim == 2) {
        fp << R"(VARIABLES = "X", "Y")";
        if (!varName.empty())
            for (const auto& it : varName) {
                fp << ", " << it;
            }
    } else {
        fp << R"(VARIABLES = "X", "Y", "Z")";
        if (!varName.empty())
            for (const auto& it : varName) {
                fp << ", " << it;
            }
    }
    fp << std::endl;
    fp << "ZONE N = " << mesh.nodeNum() << ", E = " << mesh.cellNum() << ", ";
    if (dim == 2) {
        switch (varName.size()) {
            case 0:
                fp << "VARLOCATION=([1-2]=NODAL)";
                break;
            case 1:
                fp << "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)";
                break;
            default:
                fp << "VARLOCATION=([1-2]=NODAL, [3-" << varName.size() + 2 << "]=CELLCENTERED)";
        }
        fp << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << std::endl;
    } else {
        switch (varName.size()) {
            case 0:
                fp << "VARLOCATION=([1-3]=NODAL)";
                break;
            case 1:
                fp << "VARLOCATION=([1-3]=NODAL, [4]=CELLCENTERED)";
                break;
            default:
                fp << "VARLOCATION=([1-3]=NODAL, [4-" << varName.size() + 3 << "]=CELLCENTERED)";
        }
        fp << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK" << std::endl;
    }
    fp << std::endl;
}

void writeNode(std::fstream &fp, const fvMesh& mesh) {
    const Label dim = mesh.dimension();
    const int DATA_PRECISION = 15;
    const int LINE_DATA_NUM = 30;
    int count;
    // write node
    count = 0;
    fp << std::endl << "## Node-x" << std::endl;
    for (auto &node : mesh.nodes()) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.C().x;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;


    count = 0;
    fp << std::endl << "## Node-y" << std::endl;
    for (auto &node: mesh.nodes()) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << node.C().y;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    fp << std::endl;

    if (dim == 3) {
        count = 0;
        fp << std::endl << "## Node-z" << std::endl;
        for (auto &node: mesh.nodes()) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << node.C().z;
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
        fp << std::endl;
    }
}

void writeGeom(std::fstream &fp, const fvMesh& mesh) {
    // write geom
    fp << std::endl << "## geom" << std::endl;
    for (auto &cell: mesh.cells()) {
        const int nodeNum = Geometric::nodeNum(cell.geomType());
        MESO::List<MESO::ObjectId> node(nodeNum);
        for (int i = 0; i < nodeNum; ++i) {
            node[i] = cell.nodes()[i] + 1;
        }
        switch (cell.geomType()) {
            case Geometric::Quad:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[3];
                break;
            case Geometric::Tria:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[0];
                break;
            case Geometric::Brick:
                fp << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
                   << " " << node[4] << " " << node[5] << " " << node[7] << " " << node[6];
                break;
            case Geometric::Wedge:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
                   << " " << node[3] << " " << node[4] << " " << node[5] << " " << node[5];
                break;
            case Geometric::Tetra:
                fp << " " << node[0] << " " << node[1] << " " << node[2] << " " << node[2]
                   << " " << node[3] << " " << node[3] << " " << node[3] << " " << node[3];
                break;
            case Geometric::Pyram:
                fp << " " << node[0] << " " << node[1] << " " << node[3] << " " << node[2]
                   << " " << node[4] << " " << node[4] << " " << node[4] << " " << node[4];
                break;
            default:
                logger.warn << " unsupported type." << std::endl;
                break;
        }
        fp << std::endl;
    }
}

