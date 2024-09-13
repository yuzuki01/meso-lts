#ifndef CORE_TYPE_DEF_H
#define CORE_TYPE_DEF_H

namespace MESO {
    typedef std::string String;
    typedef std::stringstream StringStream;
    typedef std::vector <String> StringList;
    typedef int ObjectId;
    typedef int ObjectType;
    typedef double Scalar;
    typedef std::string KeyString;
    typedef MESO::Math::Vector Vector;
    typedef MESO::Math::Vector Position;
    /// set
    typedef std::array<ObjectId, 2> ObjectIdSet;
    typedef std::array<Vector, 2> VectorSet;
    /// list
    typedef std::vector<ObjectId> ObjectIdList;
    typedef std::vector<ObjectIdList> GroupList;
    typedef std::vector<ObjectType> ObjectTypeList;
    typedef std::vector<Scalar> ScalarList;
    typedef std::vector<Vector> VectorList;
    /// map
    typedef std::unordered_map<String, String> StringMap;
    typedef std::unordered_map<ObjectId, ObjectId> ObjectIdMap;
    typedef std::unordered_map<String, ObjectType> ObjectTypeMap;
    typedef std::unordered_map<String, Scalar> ScalarMap;
    typedef std::unordered_map<String, Vector> VectorMap;
    typedef std::unordered_map<String, ObjectTypeList> ObjectTypeListMap;
    typedef std::unordered_map<String, ScalarList> ScalarListMap;
    typedef std::unordered_map<String, VectorList> VectorListMap;
}

#endif //CORE_TYPE_DEF_H
