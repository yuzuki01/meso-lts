#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    Time runTime(FileIO::ParamReader{"config"});
    Mesh::fvMesh mesh(
            FileIO::BasicReader("cavity.neu"),
            runTime
    );

    mesh.info();
    mesh.output();

    volScalarField a(mesh);
    // volVectorField b(mesh);

    forAll(a, ai) {
        const auto& ci = a.index()[ai];
        const auto& cell = mesh.cell(ci);
        std::cout << cell.idOnPartition() << " " << ai << std::endl;
        // a.values()[ai] = 1.0 / M_PI * exp(-(magSqr(cell.C() - Vector(0,0,0))) / 8.0);
        a.values()[ai] = mag(cell.C());
        // a[ai] = cell.V();
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
