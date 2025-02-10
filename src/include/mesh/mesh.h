#ifndef MESO_MESH_H
#define MESO_MESH_H

#include "mesh/metis.h"

#include "core/core.h"

#define namespaceMesoMesh \
    namespaceMESO \
    using namespace MESO::Mesh;


namespace MESO::Mesh {
    // pre-indication
    class Node;

    class GeomMesh;
}

#ifndef MESO_FIELD_INDICATION
#define MESO_FIELD_INDICATION
namespace MESO {
    template<typename ValueType, Label PatchType>
    class BasicField;

    enum {
        VolFlag, SurfFlag
    };

    template<typename ValueType>
    using volField = BasicField<ValueType, VolFlag>;
    using volScalarField = volField<Scalar>;
    using volVectorField = volField<Vector>;

    template<typename ValueType>
    using surfField = BasicField<ValueType, SurfFlag>;
    using surfScalarField = surfField<Scalar>;
    using surfVectorField = surfField<Vector>;

    /// LeastSquare
    class LeastSquare {
    public:
        List<Vector> dr;
        List<Scalar> weight;
        Vector Cx, Cy, Cz;

        explicit LeastSquare(const Label &neiSize)
                : dr(neiSize), weight(neiSize) {

        }
    };
}
#endif //MESO_FIELD_INDICATION

namespace MESO::Mesh {

    class BasicMeshObject {
    protected:
        const ObjectId id_;     // index
        const ObjectType GT_;   // GeomType
        const GeomMesh &mesh_;  // Owning mesh

    private:
        void operator=(const BasicMeshObject &);

    public:
        BasicMeshObject(const GeomMesh &owner,
                        const ObjectId &id,
                        const ObjectType &GT);

        ~BasicMeshObject() = default;

        /// Interfaces
        [[nodiscard]] const GeomMesh &mesh() const;

        [[nodiscard]] const ObjectId &id() const;

        [[nodiscard]] const ObjectType &geomType() const;
    };

    class NodeBasedObject : public BasicMeshObject {
    protected:
        Coordinate C_;                          // Center coordinate
        const List<ObjectId> nodes_;
        ObjectId rank_{};                       // partition
        ObjectId idOnPartition_{};              // id of partition
    private:
        void operator=(const NodeBasedObject &);

    public:
        NodeBasedObject(const GeomMesh &owner,
                        const ObjectId &id,
                        const ObjectType &GT,
                        const List<ObjectId> &nodes);

        /// Set
        void setPart(const Label &rank, const Label &idOnPart);

        /// Interfaces
        [[nodiscard]] const List<ObjectId> &nodes() const;

        [[nodiscard]] const Coordinate &C() const;

        [[nodiscard]] const Label &rank() const;

        [[nodiscard]] const Label &idOnPartition() const;
    };

    class Node : public BasicMeshObject {
    protected:
        const Coordinate C_;
    private:
        void operator=(const Node &);

    public:
        Node(const GeomMesh &owner,
             const ObjectId &id,
             const Vector &coordinate);

        ~Node() = default;

        [[nodiscard]] const Coordinate &C() const;
    };

    class Face : public NodeBasedObject {
    protected:
        Scalar S_;              // Area
        ObjectId owner_;        // owner cell
        ObjectId neighbor_;     // neighbor cell
        ObjectId patch_ = 0;    // patch id
        Vector Snv_;            // Area * normalVector
    private:
        void operator=(const Face &);

    public:
        Face(const GeomMesh &owner,
             const ObjectId &id,
             const ObjectType &geomType,
             const List<ObjectId> &nodes,
             const ObjectId &ci);

        ~Face() = default;

        void setPatch(const ObjectId &patchId);

        void setNeighbor(const Label &cellId);

        /// Interfaces
        [[nodiscard]] const Scalar &S() const;

        [[nodiscard]] const ObjectId &owner() const;

        [[nodiscard]] const ObjectId &neighbor() const;

        [[nodiscard]] const ObjectId &patch() const;

        [[nodiscard]] const Vector &Snv() const;
    };

    class Cell : public NodeBasedObject {
    protected:
        Scalar V_;                              // Volume
        List<ObjectId> faces_;                  // Faces
        List<ObjectId> neighbors_;              // Neighbor cells

    private:

        void operator=(const Cell &);

    public:
        Cell(const GeomMesh &owner,
             const ObjectId &id,
             const ObjectType &geomType,
             const List<ObjectId> &nodes
        );

        ~Cell() = default;

        void setFace(const ObjectId &idFace, const ObjectId &idOnOwner);

        void setNeighbor(const Label &nei);

        /// Interfaces
        [[nodiscard]] const Scalar &V() const;

        [[nodiscard]] const List<ObjectId> &faces() const;

        [[nodiscard]] const List<ObjectId> &neighbors() const;
    };

    class Patch {
    protected:
        const String name_;
        const GeomMesh &mesh_;
        List<ObjectId> group_;

    private:
        void operator=(const Patch &);

    public:
        Patch(const GeomMesh &owner, String name);

        ~Patch() = default;

        /// Interfaces
        void append(const ObjectId& objectId);

        void fit();

        [[nodiscard]] Label size() const;

        [[nodiscard]] const String &name() const;

        [[nodiscard]] const GeomMesh &mesh() const;

        ObjectId &operator[](const Label &index);

        const ObjectId &operator[](const Label &index) const;

        List<ObjectId> &group();

        [[nodiscard]] const List<ObjectId> &group() const;
    };

    class GeomMesh {
    protected:
        bool isParsed;
        const String name_;
        const Time &time_;

        // mesh info
        Label NNODE{};      // num of points
        Label NCELL{};      // num of cells
        Label NZONE{};      // num of zones
        Label NMARK{};      // num of marks
        Label NDFCD{};      // num of dimension
        Label NDFVL{};

        List<Node> nodes_;
        List<Cell> cells_;
        List<Patch> zones_;     // zone groups

    private:
        void operator=(const GeomMesh &);

#include "mesh/meshGeomMeshReader.h"

    public:
        explicit GeomMesh(const FileIO::BasicReader &reader,
                          const Time &time,
                          bool parse = true);

        ~GeomMesh() = default;

        /// Interfaces
        [[nodiscard]] const String &name() const;

        [[nodiscard]] const Time &time() const;

        [[nodiscard]] const Label &dimension() const;

        [[nodiscard]] const Label &nodeNum() const;

        [[nodiscard]] const Label &cellNum() const;

        [[nodiscard]] const Node &node(const ObjectId &id) const;

        [[nodiscard]] const List<Node> &nodes() const;

        [[nodiscard]] const Cell &cell(const ObjectId &id) const;

        [[nodiscard]] const List<Cell> &cells() const;
    };

    class fvMesh : public GeomMesh {
    protected:
        Label NFACE{};
        List<Face> faces_;
        List<Patch> marks_ = {
                Patch(*this, "fluid-interior")
        };
        List<Patch> partitionCellPatch_;        // cell patch
        List<Patch> partitionFacePatch_;        // face patch
        void partitionMesh(const Label &nPart);

        void generateFace();

    private:
        void operator=(const fvMesh &);

#include "mesh/meshFvMeshReader.h"

    public:
        explicit fvMesh(const FileIO::BasicReader &reader,
                        const Time &time,
                        bool parse = true);

        ~fvMesh() = default;

        /// Interfaces
        void info() const;

        [[nodiscard]] const Label &faceNum() const;

        [[nodiscard]] const Face &face(const ObjectId &id) const;

        [[nodiscard]] const List<Face> &faces() const;

        [[nodiscard]] const List<Patch> &partitionCellPatch() const;

        [[nodiscard]] const List<Patch> &partitionFacePatch() const;

        /// File IO
        void output();

        /// FVM
    protected:
        List<LeastSquare> leastSquare_;

        void initLeastSquare();

    public:
        [[nodiscard]] const List<LeastSquare> &leastSquare() const;
    };

#include "mesh/meshGeometric.h"
#include "mesh/meshBoundary.h"

}


#endif //MESO_MESH_H
