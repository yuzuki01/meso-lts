#include "mesh/mesh.h"

namespaceMesoMesh


GeomMesh::GeomMesh(const FileIO::BasicReader &reader, const Time& time, bool parse)
        : name_(reader.name()),
          time_(time),
          isParsed(parse) {
    if (parse) {
        const auto &lines = reader.read_lines();
        Label i = 0;
        auto size = static_cast<Label>(lines.size());
        parseGambitHead(i, size, lines);
        parseGambitNode(i, size, lines);
        parseGambitCell(i, size, lines);
    }
}

/**
 *  ================================================
 *  ------------------- Private --------------------
 *  ================================================
 **/


void GeomMesh::parseGambitHead(Label &i, const Label &size,
                               const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        data = Utils::split(lines[i]);
        if (data[0] == "NUMNP") {
            data = Utils::split(lines[++i]);
            NNODE = std::stoi(data[0]);
            NCELL = std::stoi(data[1]);
            NZONE = std::stoi(data[2]);
            zones_.reserve(NZONE);
            NMARK = std::stoi(data[3]) + 1;
            NDFCD = std::stoi(data[4]);
            NDFVL = std::stoi(data[5]);
        }
        if (data[0] == "ENDOFSECTION") break;
    }
    ++i;
    logger.debug << "Read mesh head - Done" << std::endl;
}

void GeomMesh::parseGambitNode(Label &i, const Label &size,
                               const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        data = Utils::split(lines[i]);
        if (data[0] == "NODAL" and data[1] == "COORDINATES") {
            ++i;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                // node
                Label ni = std::stoi(data[0]) - 1;
                Scalar nx = 0, ny = 0, nz = 0;
                if (data.size() >= 3) {
                    nx = std::stod(data[1]);
                    ny = std::stod(data[2]);
                }
                if (data.size() >= 4) {
                    nz = std::stod(data[3]);
                }
                nodes_.emplace_back(*this, ni,
                                    Vector(nx, ny, nz));
                ++i;
            }
        }
        if (data[0] == "ENDOFSECTION") break;
    }
    ++i;
    logger.debug << "Read mesh nodes - Done" << std::endl;
}

void GeomMesh::parseGambitCell(Label &i, const Label &size,
                               const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        data = Utils::split(lines[i]);
        if (data[0] == "ELEMENTS/CELLS") {
            ++i;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                // cell
                ObjectId ci = std::stoi(data[0]) - 1;
                ObjectType geomType = std::stoi(data[1]);
                Label nNum = std::stoi(data[2]);
                auto dataLen = static_cast<Label>(data.size());
                List<ObjectId> node_list;
                for (Label j = 3; j < dataLen; ++j) {
                    ObjectId ni = std::stoi(data[j]) - 1;
                    node_list.push_back(ni);
                }
                if (node_list.size() < nNum) {
                    ++i;
                    data = Utils::split(lines[i]);
                    for (const auto &it: data) {
                        ObjectId ni = std::stoi(it) - 1;
                        node_list.push_back(ni);
                    }
                }
                cells_.emplace_back(*this, ci,
                                    geomType, node_list);
                ++i;
            }
            if (data[0] == "ENDOFSECTION") break;
        }
    }
    ++i;
    logger.debug << "Read mesh cells - Done" << std::endl;
}


/**
 *  ================================================
 *  ------------------ Interfaces ------------------
 *  ================================================
 **/


const String &GeomMesh::name() const {
    return name_;
}

const Time &GeomMesh::time() const {
    return time_;
}

const Label &GeomMesh::dimension() const {
    return NDFCD;
}

const Label &GeomMesh::nodeNum() const {
    return NNODE;
}

const Label &GeomMesh::cellNum() const {
    return NCELL;
}

const Node &GeomMesh::node(const ObjectId &id) const {
    return nodes_[id];
}

const List<Node> &GeomMesh::nodes() const {
    return nodes_;
}

const Cell &GeomMesh::cell(const ObjectId &id) const {
    return cells_[id];
}

const List<Cell> &GeomMesh::cells() const {
    return cells_;
}
