#include "mesh/mesh.h"

namespaceMesoMesh

using namespace Geometric;

Map<List<ObjectId>> fgMap = {
        {Quad,  {Edge, Edge, Edge, Edge}},
        {Tria,  {Edge, Edge, Edge}},
        {Brick, {Quad, Quad, Quad, Quad, Quad, Quad}},
        {Wedge, {Quad, Quad, Quad, Tria, Tria}},
        {Tetra, {Tria, Tria, Tria, Tria}},
        {Pyram, {Quad, Tria, Tria, Tria, Tria}}
};

Map<List<List<ObjectId>>> fniMap = {
        {Quad,  {{0, 1},       {1, 2},       {2, 3},       {3, 0}}},
        {Tria,  {{0, 1},       {1, 2},       {2, 0}}},
        {Brick, {{0, 1, 5, 4}, {1, 3, 7, 5}, {3, 2, 6, 7}, {2, 0, 4, 6}, {1, 0, 2, 3}, {4, 5, 7, 6,}}},
        {Wedge, {{0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {0, 2, 1},    {3, 4, 5}}},
        {Tetra, {{1, 0, 2},    {0, 1, 3},    {1, 2, 3},    {2, 0, 3}}},
        {Pyram, {{0, 2, 3, 1}, {0, 1, 4},    {1, 3, 4},    {3, 2, 4},    {2, 0, 4}}}
};

KeyString Geometric::generateKey(const List<ObjectId> &nodes) {
    const int len = int(nodes.size());
    /// find where is min
    int pos = 0;
    int min = nodes[pos];
    for (int i = 0; i < len; ++i) {
        if (nodes[i] < min) {
            pos = i;
            min = nodes[i];
        }
    }
    /// direction
    bool direction;
    int left_pos, right_pos;
    left_pos = (pos == 0) ? (len - 1) : (pos - 1);
    right_pos = (pos == len - 1) ? 0 : (pos + 1);
    direction = nodes[right_pos] < nodes[left_pos];
    /// append key
    int count = 0;
    std::stringstream ss;
    while (count < len) {
        if (count == 0) ss << nodes[pos]; else ss << ">" << nodes[pos];
        pos = direction ? (pos + 1) : (pos - 1);
        /// make valid range
        if (pos < 0) pos = len - 1;
        if (pos > len - 1) pos = 0;
        count++;
    }
    return ss.str();
}


void fvMesh::generateFace() {
    Dict<ObjectId> fmEdge, fmQuad, fmTria;
    for (auto &cell: cells_) {
        for (int idOnOwner = 0; idOnOwner < Geometric::faceNum(cell.geomType()); ++idOnOwner) {
            /// range faces on cell
            const auto &geomType = fgMap[cell.geomType()][idOnOwner];
            List<ObjectId> nodeList;
            for (const auto &it: fniMap[cell.geomType()][idOnOwner]) {
                nodeList.push_back(cell.nodes()[it]);
            }
            String faceKey = Geometric::generateKey(nodeList);
            Dict<ObjectId> *fmp;
            switch (geomType) {
                case Geometric::Edge:
                    fmp = &fmEdge;
                    break;
                case Geometric::Quad:
                    fmp = &fmQuad;
                    break;
                case Geometric::Tria:
                    fmp = &fmTria;
                    break;
                default:
                    logger.error << "MESO::fvMesh::generateFace() caught unexpected value" << std::endl;
                    FATAL_ERROR_THROW;
            }
            Dict<ObjectId> &fm = *fmp;
            /// check face
            auto it = fm.find(faceKey);
            if (it == fm.end()) {
                /// not exist
                auto faceId = Label(faces_.size());
                faces_.emplace_back(*this, faceId,
                                    geomType, nodeList,
                                    cell.id());
                fm[faceKey] = faceId;
                auto &face = faces_.back();
                cell.setFace(face.id(), idOnOwner);
            } else {
                /// exist
                auto &face = faces_[it->second];
                face.setNeighbor(cell.id());
                cell.setFace(face.id(), idOnOwner);
            }
        }
    }
    NFACE = Label(faces_.size());
    /// Construct Cell-neighbor
    for (const auto& face : faces_) {
        if (face.owner() == face.neighbor()) continue;
        auto &owner = cells_[face.owner()];
        auto &neighbor = cells_[face.neighbor()];
        owner.setNeighbor(face.neighbor());
        neighbor.setNeighbor(face.owner());
    }
}
