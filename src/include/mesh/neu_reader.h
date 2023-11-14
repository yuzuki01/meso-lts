/**
 * included by mesh.h
 */

class neuReader : public BasicReader {
protected:
    const std::string end_mark = "ENDOFSECTION";
public:
    explicit neuReader(const std::string &file_path) : BasicReader("neuReader", file_path) {};

    TP_key void parse(MESH::Mesh<key_type> &mesh);
};
