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

        [[nodiscard]] const GeomMesh &mesh() const;

        [[nodiscard]] const ObjectId &id() const;

        [[nodiscard]] const ObjectType &geomType() const;
    };

    class NodeBasedObject : public BasicMeshObject {
    protected:
        Vector C_;              // Center coordinate
        const List<ObjectId> nodes_;
    private:
        void operator=(const NodeBasedObject &);

    public:
        NodeBasedObject(const GeomMesh &owner,
                        const ObjectId &id,
                        const ObjectType &GT,
                        const List<ObjectId> &nodes);

        [[nodiscard]] const List<ObjectId> &nodes() const;

        [[nodiscard]] const Vector &C() const;
    };

    class Node : public BasicMeshObject {
    protected:
        const Vector C_;    // coordinate
    private:
        void operator=(const Node &);

    public:
        Node(const GeomMesh &owner,
             const ObjectId &id,
             const Vector &coordinate);

        ~Node() = default;

        [[nodiscard]] const Vector &C() const;
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

        void setNeighbor(const Label& cellId);

        /// Interfaces

        [[nodiscard]] const Scalar &S() const;

        [[nodiscard]] const ObjectId &owner() const;

        [[nodiscard]] const ObjectId &neighbor() const;

        [[nodiscard]] const ObjectId &patch() const;

        [[nodiscard]] const Vector &Snv() const;
    };

    class Cell : public NodeBasedObject {
    protected:
        Scalar V_;                  // Volume
        List<ObjectId> faces_;      // Faces

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

        /// Interfaces

        [[nodiscard]] const Scalar &V() const;

        [[nodiscard]] const List<ObjectId> &faces() const;
    };

    class Patch {
    protected:
        const String name;
        const GeomMesh &mesh_;
        List<ObjectId> group_;

    private:
        void operator=(const Patch &);

    public:
        Patch(const GeomMesh &owner, String name);

        ~Patch() = default;

        [[nodiscard]] const GeomMesh &mesh() const;

        List<ObjectId> &group();

        [[nodiscard]] const List<ObjectId> &group() const;
    };

    class GeomMesh {
    protected:
        bool isParsed;
        const String name_;
        const Time time_;

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
        List<Patch> parts_;     // cell partitions
        void partition(const Label &nPart);

    private:
        void operator=(const fvMesh &);

#include "mesh/meshFvMeshReader.h"

        void generateFace();

    public:
        explicit fvMesh(const FileIO::BasicReader &reader,
                        const Time &time,
                        bool parse = true);

        ~fvMesh() = default;

        /// Interfaces
        [[nodiscard]] const Label &faceNum() const;

        [[nodiscard]] const Face &face(const ObjectId &id) const;

        [[nodiscard]] const List<Face> &faces() const;

        /// File IO
        void output();
    };

#include "mesh/meshGeometric.h"

}

#endif //MESO_MESH_H
