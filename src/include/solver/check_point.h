template <class Solver>
class CheckPoint {
private:
    Solver &solver;

public:
    explicit CheckPoint(Solver &_solver) : solver(_solver) {};
    void init_field(const PhysicalVar::MacroVars &_var);
    void init_from_file(const std::string &file_path);
    void write_to_file(const std::string &file_path);
};

static PhysicalVar::MacroVars strvec_to_macro_vars(const string_vector &str_vec);
