//
// Created by Administrator on 2024/5/29.
//

#ifndef SOLVER_CONFIG_H
#define SOLVER_CONFIG_H


namespace MESO::Solver {

    extern std::unordered_map<std::string, int> mark_type_map;

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
    }

    class Patch {
    private:
        ObjectTypeMap type;
        ObjectTypeMap ints;
        ScalarMap scalars;
        VectorMap vectors;
    public:
        void set_value(StringList &string_list);

        ObjectType get_type(const String &key);

        ObjectType get_int(const String &key);

        Scalar get_scalar(const String &key);

        Vector get_vector(const String &key);
    };

    struct Mark {
        std::string name;
        ObjectId id;
        ObjectType type;
        std::string type_name;
        Patch patch;
        /*
        double density;             // inlet density
        double temperature=0.0;         // wall temperature, flow temperature
        double pressure;            // total pressure
        Vector velocity;            // inlet velocity, wall velocity
        double slip_wall_alpha=0.0;
        int direction=0;
         */
    };

    struct Group {
        std::string name;
        ObjectId id;
        Patch patch;
        /*
        double density;
        double temperature;
        Vector velocity;
         */
    };

    typedef std::string ConfigKey;

    class Config {
    public:
        std::unordered_map<std::string, std::string> settings;
        std::unordered_map<std::string, Group> groups;
        std::unordered_map<std::string, Mark> marks;

        Config() = default;

        explicit Config(const std::string &file_path);

        /**
         * T = [int, double, bool, std::string]
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
