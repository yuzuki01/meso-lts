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
    template<class T>
    using Set = std::array<T, 2>;
    /// list
    template<class T>
    using List = std::vector<T>;
    /// Map
    template<class T>
    using Map = std::unordered_map<ObjectId, T>;
    /// Dict
    template<class T>
    using Dict = std::unordered_map<String, T>;
}

#endif //CORE_TYPE_DEF_H
