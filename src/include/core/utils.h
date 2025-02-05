#ifndef MESO_UTILS_H
#define MESO_UTILS_H

namespace MESO::Utils {
    StringList split(const std::string &_str);

    void print_names_and_values(const StringList &names, const List<Scalar> &values);

    int mkdir(const String &dir_name);

    bool is_converged(const List<Scalar> &residual_list, Scalar limit);

    template<class T>
    void output_list(const String &file_path, List<T> &data);

    template<class T>
    List<T> read_np_file(const String &file_path);

    /// receive python3 script prints
    StringList exec_script(const String &script);
}

#endif //MESO_UTILS_H
