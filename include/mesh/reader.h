#include <utility>

#ifndef MESH_READER_H
#define MESH_READER_H

namespace MESO::fvmMesh::Reader {
    class Gambit;
}

class MESO::fvmMesh::Reader::Gambit : public MESO::FileReader::BasicReader {
public:
    explicit Gambit(MESO::String file) : MESO::FileReader::BasicReader(std::move(file)) {};
    Mesh parse_mesh_strings();
};

#endif //MESH_READER_H
