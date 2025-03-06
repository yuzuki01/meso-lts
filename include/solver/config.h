//
// Created by Administrator on 2024/5/29.
//

#ifndef SOLVER_CONFIG_H
#define SOLVER_CONFIG_H


namespace MESO::Solver {

    extern const Dict<int> mark_type_map;

    enum BoundaryType {
        fluid_interior,
        farfield_inlet,
        pressure_inlet,
        freestream_inlet,
        farfield_outlet,
        pressure_outlet,
        wall,
        slip_wall,
        symmetry,
    };

    namespace PatchType {
        const ObjectType zeroGradient = 0;
        const String zeroGradient_str = "zeroGradient";
        const ObjectType fixedValue = 1;
        const String fixedValue_str = "fixedValue";
        const ObjectType calculated = 2;
        const String calculated_str = "calculated";
        const ObjectType fromFile = 3;
        const String fromFile_str = "fromFile";
        const ObjectType switchType = 4;
        const String switchType_str = "switch";
    }

    class Patch {
    private:
        Dict<ObjectType> type;
        Dict<ObjectType> ints;
        Dict<bool> switches = {
                /// default value
                {"mask", false}
        };
        Dict<Scalar> scalars;
        Dict<Vector> vectors;
        Dict<List<ObjectType>> ints_list;
        Dict<List<Scalar>> scalars_list;
        Dict<List<Vector>> vectors_list;
    public:
        void set_value(StringList &string_list);

        ObjectType get_type(const String &key);

        bool get_switch(const String &key);

        ObjectType get_int(const String &key);

        Scalar get_scalar(const String &key);

        Vector get_vector(const String &key);

        ObjectType get_file_int(const String &key, ObjectId id);

        Scalar get_file_scalar(const String &key, ObjectId id);

        Vector get_file_vector(const String &key, ObjectId id);
    };

    struct Mark {
        String name;
        ObjectId id;
        ObjectType type;
        String type_name;
        Patch patch;
    };

    struct Group {
        String name;
        ObjectId id;
        Patch patch;
    };

    typedef String ConfigKey;

    class Config {
    private:
        const String file_path;
    public:
        Dict<String> settings;
        Dict<Group> groups;
        Dict<Mark> marks;

        Config() = default;

        explicit Config(const String &file_path);

        void update_config(bool update_patch=true);

        /**
         * T = [int, double, bool, String]
         * 用法:
         *  auto is_limiter_open = config.get<bool>("limiter-switch", true);
         **/
        template<class T>
        T get(const ConfigKey &key, T default_value, bool print_error = true);

        Group &get_cell_group(MESO::fvmMesh::Cell &cell, MESO::fvmMesh::Mesh &mesh);

        Mark &get_face_group(MESO::fvmMesh::Face &face, MESO::fvmMesh::Mesh &mesh);

        void info();
    };
}


#endif //SOLVER_CONFIG_H
