/**
 * included by core.h
 */

/// 文本读取基类
class BasicReader {
protected:
    const int print_line_num = 20;
    bool file_open = false;
    std::string name{}, path{};
    std::vector<std::string> lines{};
    Logger logger{};
public:
    explicit BasicReader(const std::string &reader_name, const std::string &file_path);

    BasicReader() = default;

    ~BasicReader() = default;

    bool is_file_open() const;

    void print();
};
