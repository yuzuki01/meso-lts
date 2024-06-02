#include <utility>

#ifndef MESH_READER_H
#define MESH_READER_H

namespace MESO::Mesh::Reader {
    class BasicReader;
    class Gambit;
}

class MESO::Mesh::Reader::BasicReader {
public:
    String file;
    explicit BasicReader(String file) : file(std::move(file)) {};

    [[nodiscard]] StringList read_lines() const;

    virtual Mesh read() { throw std::invalid_argument("<BasicReader> read() has not been override."); };
};

class MESO::Mesh::Reader::Gambit : public MESO::Mesh::Reader::BasicReader {
public:
    explicit Gambit(MESO::String file) : MESO::Mesh::Reader::BasicReader(std::move(file)) {};
    Mesh read() override;
};

#endif //MESH_READER_H
