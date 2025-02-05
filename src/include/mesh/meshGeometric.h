/**
 * File: mesh.h
 *  namespace MESO::fvMesh {
 *      #include "mesh/meshGeometric.h"
 *  }
 **/


#ifndef MESO_MESHGEOMTRIC_H
#define MESO_MESHGEOMTRIC_H

namespace Geometric {
    // dim-0
    const ObjectType Node = 0;
    // dim-1
    const ObjectType Edge = 1;
    // dim-2
    const ObjectType Quad = 2;
    const ObjectType Tria = 3;
    // dim-3
    const ObjectType Brick = 4;
    const ObjectType Wedge = 5;
    const ObjectType Tetra = 6;
    const ObjectType Pyram = 7;

    Label nodeNum(const ObjectType& geomType);
    Label faceNum(const ObjectType& geomType);
    Vector getCoordinate(const GeomMesh &mesh,
                         const List<ObjectId> &nodes);
    Scalar getArea(const GeomMesh &mesh,
                   const List<ObjectId>& nodes,
                   const ObjectType& GT);
    Vector getSurfaceVector(const GeomMesh &mesh,
                            const List<ObjectId>& nodes,
                            const ObjectType& GT,
                            const ObjectId &owner,
                            const Vector& faceC,
                            const Scalar& Sf);
    Scalar getVolume(const GeomMesh &mesh,
                     const List<ObjectId>& nodes,
                     const ObjectType &GT);
    // fvMeshGenerateFace
    KeyString generateKey(const List<ObjectId>& nodes);
}

#endif //MESO_MESHGEOMTRIC_H
