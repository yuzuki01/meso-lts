template <class Solver>
class CheckPoint {
private:
    Solver &solver;
    const std::string end_mark = ":end";
    const std::string cell_mark = "cell:";
public:
    explicit CheckPoint(Solver &_solver) : solver(_solver) {};
    void init_field(const Physical::MacroVars &_var);
    void init_from_file(const std::string &file_path);
    void write_to_file(const std::string &file_path);
};
