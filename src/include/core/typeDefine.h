#ifndef CORE_TYPE_DEF_H
#define CORE_TYPE_DEF_H

namespace MESO {
    namespace Math {
        class Vector;
    }

    typedef std::string String;
    typedef std::stringstream StringStream;
    typedef std::vector <String> StringList;
    typedef int Label;
    typedef int ObjectId;
    typedef int ObjectType;
    typedef double Scalar;
    typedef std::string KeyString;
    typedef MESO::Math::Vector Vector;
    typedef MESO::Math::Vector Coordinate;
    /// set
    template<class T>
    using Set = std::array<T, 2>;
    /// list
    template<class T>
    using List = std::vector<T>;
    /// Map
    template<class T>
    using Map = std::unordered_map<Label, T>;
    /// Dict
    template<class T>
    using Dict = std::unordered_map<String, T>;
    /// Enum
    enum {X_COMPONENT, Y_COMPONENT, Z_COMPONENT};
}

#endif //CORE_TYPE_DEF_H
