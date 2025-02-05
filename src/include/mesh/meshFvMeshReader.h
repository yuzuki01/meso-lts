#include "mesh/meshGeomMeshReader.h"
private:
    void parseGambitZone(Label &i, const Label &size,
                         const StringList &lines);
    void parseGambitMark(Label &i, const Label &size,
                         const StringList &lines);
