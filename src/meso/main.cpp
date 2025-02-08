#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    Time runTime(FileIO::ParamReader{"config"});
    Mesh::fvMesh mesh(
            FileIO::BasicReader("dvs.neu"),
            runTime
    );

    mesh.info();
    mesh.output();

    volScalarField a(mesh);
    // volVectorField b(mesh);

    forAll(a, ai) {
        const auto& ci = a.index()[ai];
        const auto& cell = mesh.cell(ci);
        // a[ai] = 1.0 / M_PI * exp(-(magSqr(cell.C() - Vector(0.25,0,0))));
        a[ai] = cell.V();
        // b[ai] = cell.C();
    }

    a.output("a");
    // b.output("b");

    auto grad = fvm::grad(a);

    grad.output("grad");

    MPI::Finalize();

    /// run with no params
    logger.warn << "type \"./meso-mpi -help\" to get help." << std::endl;
    return 0;
}
