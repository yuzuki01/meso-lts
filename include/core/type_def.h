#ifndef CORE_TYPE_DEF_H
#define CORE_TYPE_DEF_H

namespace MESO {
    typedef std::string String;
    typedef std::vector <String> StringList;
    typedef int ObjectId;
    typedef int ObjectType;
    typedef double Scalar;
    typedef std::string KeyString;
    typedef MESO::Math::Vector Vector;
    typedef MESO::Math::Matrix Matrix;
    typedef MESO::Math::Vector Position;
    typedef std::array<ObjectId, 2> ObjectIdSet;
    typedef std::array<Vector, 2> VectorSet;
    typedef std::vector<ObjectId> ObjectIdList;
    typedef std::vector<ObjectIdList> GroupList;
    typedef std::vector<Scalar> ScalarList;
    typedef std::vector<Vector> VectorList;
}

#endif //CORE_TYPE_DEF_H
