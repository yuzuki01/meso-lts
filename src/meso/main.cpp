#include "meso.h"

namespaceMESO


int main(int argc, char **argv) {
    ArgParser parser(argc, argv);

    if (parser.parse_switch("debug")) {
        logger.level = -1;
    }
    logger.debug << "Running in debug mode." << std::endl;

    MPI::Initialize(&argc, &argv);

    Time runTime(FileIO::ParamReader {"config"});
    Mesh::fvMesh mesh(
            FileIO::BasicReader("cavity.neu"),
            runTime
            );
    logger.info << "nFace=" << mesh.faces().size() << std::endl;

    mesh.output();

    MPI::Finalize();

    /// run with no params
    logger.warn << "type \"./meso-mpi -h\" to get help." << std::endl;
    return 0;
}
