#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    MPI::Initialize(&argc, &argv);

    Time runTime(FileIO::ParamReader{"config"});
    Mesh::fvMesh mesh(
            FileIO::BasicReader("cavity.neu"),
            runTime
    );
    logger.info << "nFace=" << mesh.faces().size() << std::endl;

    mesh.info();
    mesh.output();

    volVectorField a(mesh);
    volScalarField b(mesh);

    MPI::Finalize();

    /// run with no params
    logger.warn << "type \"./meso-mpi -help\" to get help." << std::endl;
    return 0;
}
