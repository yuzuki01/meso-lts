#ifndef MESO_DVS_READER_H
#define MESO_DVS_READER_H

namespace MESO::DVS {

    struct GaussHermiteParams {
        int dimension;
        double RT;

        GaussHermiteParams(int dimension, double RT) : dimension(dimension), RT(RT) {};
    };

    struct HalfRangeGHParams {
        int dimension;
        double RT;

        HalfRangeGHParams(int dimension, double RT) : dimension(dimension), RT(RT) {};
    };

    struct NewtonCotesParams {
        int dimension;

        explicit NewtonCotesParams(int dimension) : dimension(dimension) {};
    };

    template<typename T>
    fvmMesh::Mesh generate_dvs(const String &file_path, T &params);

    /// Gauss-Hermite
    template<>
    fvmMesh::Mesh generate_dvs(const String &file_path, GaussHermiteParams &params);

    /// half-range Gauss-Hermite
    template<>
    fvmMesh::Mesh generate_dvs(const String &file_path, HalfRangeGHParams &params);

    /// Newton-Cotes
    template<>
    fvmMesh::Mesh generate_dvs(const String &file_path, NewtonCotesParams &params);

    /// Sparse
    fvmMesh::Mesh read_mesh_file(const String &file_path);
}

#endif //MESO_DVS_READER_H
