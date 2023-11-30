/**
 * included by core.h
 */

void meso_init();

string_vector split(const std::string &_str);

bool str_vec_cmp(const string_vector &_str_vec, const std::string &reference);


/// println 全局函数
void info_println(const std::string &context);

void note_println(const std::string &context);

void warn_println(const std::string &context);

void error_println(const std::string &context);

void debug_println(const std::string &context);

void green_println(const std::string &context);

void highlight_println(const std::string &context);

void data_int_println(const std::vector<std::string> &names, const std::vector<int> &values, int width=15);

void data_double_println(const std::vector<std::string> &names, const std::vector<double> &values, int width=15);

void data_sci_double_println(const std::vector<std::string> &names, const std::vector<double> &values, int width=15);


/// file
bool create_dir(const std::string &path);
