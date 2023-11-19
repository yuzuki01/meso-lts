/**
 * included by mesh.h
 */

class su2Reader : public BasicReader {
protected:
    const std::string ndime_mark = "NDIME=",
        npoin_mark = "NPOIN=", nelem_mark = "NELEM=",
        nmark_mark = "NMARK=", markertag_mark = "MARKER_TAG=", markerelems_mark = "MARKER_ELEMS=";
public:
    explicit su2Reader(const std::string &file_path) : BasicReader("su2Reader", file_path) {};

    TP_mesh void parse(mesh_type &mesh);
};
