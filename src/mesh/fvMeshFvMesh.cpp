#include "mesh/mesh.h"

namespaceMesoMesh

fvMesh::fvMesh(const String &filePath, const Time &time, bool parse)
        : GeomMesh(filePath, time, false) {
    isParsed = parse;
    if (parse) {
        FileIO::BasicReader reader(filePath);
        const auto &lines = reader.read_lines();
        Label i = 0;
        auto size = static_cast<Label>(lines.size());
        parseGambitHead(i, size, lines);
        parseGambitNode(i, size, lines);
        parseGambitCell(i, size, lines);
        parseGambitZone(i, size, lines);
        parseGambitMark(i, size, lines);
        // partition
        partitionMesh(MPI::processorNum);
        // init LeastSquare here
        initLeastSquare();
    }
}

/**
 *  ================================================
 *  ------------------- Private --------------------
 *  ================================================
 **/


void fvMesh::parseGambitHead(Label &i, const Label &size,
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

void fvMesh::parseGambitNode(Label &i, const Label &size,
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

void fvMesh::parseGambitCell(Label &i, const Label &size,
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
    /// generate face here
    generateFace();
    logger.debug << "Read mesh cells - Done" << std::endl;
}

void fvMesh::parseGambitZone(Label &i, const Label &size,
                             const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        data = Utils::split(lines[i]);
        if (data[0] == "ELEMENT" and data[1] == "GROUP") {
            ++i;
            data = Utils::split(lines[i]);
            int group_id = std::stoi(data[1]) - 1;
            ++i;
            data = Utils::split(lines[i]);
            std::string group_name = data[0];
            zones_.emplace_back(*this, group_name);
            i += 2;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                /// parse
                for (const auto &it: data) {
                    int cell_id = std::stoi(it) - 1;
                    zones_[group_id].group().push_back(cell_id);
                }
                ++i;
            }
        }
        if (data[0] == "BOUNDARY" and data[1] == "CONDITIONS") {
            --i;
            break;
        }
    }
    ++i;
    logger.debug << "Read mesh zones - Done" << std::endl;
}

void fvMesh::parseGambitMark(Label &i, const Label &size,
                             const StringList &lines) {
    StringList data;
    for (; i < size; ++i) {
        data = Utils::split(lines[i]);
        if (data[0] == "BOUNDARY" and data[1] == "CONDITIONS") {
            ++i;
            data = Utils::split(lines[i]);
            String bc_name = data[0];
            Label nFace = std::stoi(data[2]);
            auto idMark = static_cast<Label>(marks_.size());
            marks_.emplace_back(*this, bc_name);
            auto &mark_ = marks_.back();
            mark_.group().reserve(nFace);
            ++i;
            while (i < size) {
                data = Utils::split(lines[i]);
                if (data[0] == "ENDOFSECTION") break;
                ObjectId idCell = std::stoi(data[0]) - 1;
                ObjectId idFaceOfOwner = std::stoi(data[2]) - 1;
                auto &owner = cells_[idCell];
                auto &face_ = faces_[owner.faces()[idFaceOfOwner]];
                face_.setPatch(idMark);
                ++i;
            }
        }
    }
    for (auto &face_: faces_) {
        marks_[face_.patch()].group().push_back(face_.id());
    }
    logger.debug << "Read mesh marks - Done" << std::endl;
}


/**
 *  ================================================
 *  ------------------ Protected -------------------
 *  ================================================
 **/

void fvMesh::partitionMesh(const Label &nPart) {
    if (nPart == 1) {
        partitionCellPatch_.reserve(1);
        partitionCellPatch_.emplace_back(*this, "cellPartition");
        auto &patchCell = partitionCellPatch_.back();
        partitionFacePatch_.reserve(1);
        patchCell.group().reserve(NCELL);
        partitionFacePatch_.emplace_back(*this, "facePartition");
        auto &patchFace = partitionFacePatch_.back();
        patchFace.group().reserve(NFACE);
        for (auto &cell_: cells_) {
            patchCell.append(cell_.id());
            cell_.setPart(0, cell_.id());
        }
        for (auto &face_: faces_) {
            patchFace.append(face_.id());
            face_.setPart(0, face_.id());
        }
        logger.debug << "No need to partitionMesh mesh " << name_ << std::endl;
        return;
    }
    // call METIS
    auto numCells = static_cast<idx_t>(NCELL);
    auto numParts = static_cast<idx_t>(nPart);
    typedef std::vector<idx_t> MetisIdList;
    MetisIdList partitionCells(numCells);

    if (MPI::isMainProcessor()) {
        List<List<ObjectId>> neighbors(numCells);
        for (int ci = 0; ci < numCells; ++ci) {
            neighbors[ci] = cells_[ci].neighbors();
        }

        MetisIdList xadj(numCells + 1);
        MetisIdList adjncy;

        xadj[0] = 0;
        for (int i = 0; i < numCells; ++i) {
            Cell &cell = cells_[i];
            xadj[i + 1] = xadj[i] + static_cast<Label>(neighbors[cell.id()].size());
            for (auto neighbor: neighbors[cell.id()]) {
                adjncy.push_back(static_cast<idx_t>(neighbor));
            }
        }

        idx_t objval;
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        METIS_PartGraphKway(&numCells, &numParts, xadj.data(), adjncy.data(),
                            nullptr, nullptr, nullptr, &numParts, nullptr,
                            nullptr, options, &objval, partitionCells.data());
    }

    MPI::Bcast(partitionCells, MPI::mainRank
    );

    MPI::Barrier();
    // METIS to MESO
    partitionCellPatch_.reserve(nPart);
    partitionFacePatch_.reserve(nPart);
    for (int i = 0; i < nPart; ++i) {
        StringStream partName;
        partName << "Part:" << i;
        partitionCellPatch_.emplace_back(*this, partName.str());
        partitionFacePatch_.emplace_back(*this, partName.str());
    }
    for (int i = 0; i < numCells; ++i) {
        auto &partPatch = partitionCellPatch_[partitionCells[i]];
        cells_[i].setPart(partitionCells[i],
                          partPatch.size());
        partPatch.append(i);
    }
    for (auto &part: partitionCellPatch_) {
        part.fit();
    }
    for (int i = 0; i < NFACE; ++i) {
        auto &face = faces_[i];
        const auto &owner = cells_[face.owner()];
        auto &partPatch = partitionFacePatch_[owner.rank()];
        face.setPart(owner.rank(),
                     partPatch.size());
        partPatch.append(face.id());
    }
    for (auto &part: partitionFacePatch_) {
        part.fit();
    }
    // Reset Interior Faces
    auto& interiorPatch = marks_[0];
    Patch processorPatch(*this, "processor");
    List<ObjectId> interiorFaces;
    for (const Label &fi : interiorPatch.group()) {
        const auto& face = faces_[fi];
        const auto& own = cells_[face.owner()];
        const auto& nei = cells_[face.neighbor()];
        if (partitionCells[own.id()] == partitionCells[nei.id()]) {
            interiorFaces.push_back(fi);
        } else {
            processorPatch.group().push_back(fi);
        }
    }
    if (processorPatch.size() > 0) {
        interiorPatch.group() = interiorFaces;
        marks_.push_back(processorPatch);
        marks_.shrink_to_fit();
    }
    for (auto &part: marks_) {
        part.fit();
    }
    logger.debug << "Mesh {" << name_ << "} was partitioned into " <<
                 Label(partitionCellPatch_.size()) << " part(s)" <<
                 std::endl;
}


/**
 *  ================================================
 *  ------------------ Interfaces ------------------
 *  ================================================
 **/

void fvMesh::info() const {
    if (not MPI::isMainProcessor()) return;
    logger.info << "\nfvMesh:"
                   "\n    name: " << name_ <<
                   "\n    num of nodes: " << NNODE <<
                   "\n    num of cells: " << NCELL <<
                   "\n    num of faces: " << NFACE <<
                   "\n    num of marks: " << NMARK <<
                   std::endl;
    // zone info
    if (not zones_.empty()) {
        logger.info << "  fvMesh::Zones:";
        for (int i = 0; i < zones_.size(); ++i) {
            const auto& zonePatch = zones_[i];
            logger.info << "\n    " << i << " - " << zonePatch.name()
                        << "\n        num of faces: " << zonePatch.size();
        }
        logger.info << std::endl;
    }
    // mark info
    if (not marks_.empty()) {
        logger.info << "  fvMesh::Marks:";
        for (int i = 0; i < marks_.size(); ++i) {
            const auto& markPatch = marks_[i];
            logger.info << "\n    " << i << " - " << markPatch.name()
                        << "\n        num of faces: " << markPatch.size();
        }
        logger.info << std::endl;
    }
    logger.info << std::endl;
}

const Label &fvMesh::faceNum() const {
    return NFACE;
}

const Face &fvMesh::face(const ObjectId &id) const {
    return faces_[id];
}

const List<Face> &fvMesh::faces() const {
    return faces_;
}

const Patch &fvMesh::mark(const ObjectId &id) const {
    return marks_[id];
}

const List<Patch> &fvMesh::marks() const {
    return marks_;
}

const List<Patch> &fvMesh::partitionCellPatch() const {
    return partitionCellPatch_;
}

const List<Patch> &fvMesh::partitionFacePatch() const {
    return partitionFacePatch_;
}
