//
// Created by Administrator on 2024/5/29.
//

#ifndef SOLVER_CONFIG_H
#define SOLVER_CONFIG_H


namespace MESO::Solver {

    extern std::unordered_map<std::string, int> mark_type_map;

    enum {fluid_interior, inlet, outlet, wall};

    struct Mark {
        std::string name;
        ObjectId id;
        ObjectType type;
        std::string type_name;
        double density;
        double temperature;
        Vector velocity;
    };

    struct Group {
        std::string name;
        ObjectId id;
        double density;
        double temperature;
        Vector velocity;
    };

    typedef std::string ConfigKey;

    class Config {
    public:
        std::unordered_map<std::string, std::string> settings;
        std::unordered_map<std::string, Group> groups;
        std::unordered_map<std::string, Mark> marks;

        Config() = default;

        explicit Config(const std::string &file_path);

        template<class T>
        /**
         * T = [int, double, bool, std::string]
         * 用法:
         *  auto is_limiter_open = config.get<bool>("limiter-switch", true);
         **/
        T get(const ConfigKey &key, T default_value, bool print_error = true);

        Group & get_cell_group(MESO::Mesh::Cell &cell, MESO::Mesh::Zone &mesh);
        Mark & get_face_group(MESO::Mesh::Face &face, MESO::Mesh::Zone &mesh);

        void info();
    };
}


#endif //SOLVER_CONFIG_H
