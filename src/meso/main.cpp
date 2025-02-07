#include "meso.h"


int main(int argc, char **argv) {

#include "mesoInit.h"

#include "mesoHelp.h"

    MPI::Initialize(&argc, &argv);

    Time runTime(FileIO::ParamReader{"config"});
    Mesh::fvMesh mesh(
            FileIO::BasicReader("dvs.neu"),
            runTime
    );

    mesh.info();
    mesh.output();

    volScalarField a(mesh);

    forAll(a, ai) {
        const auto& ci = a.index()[ai];
        const auto& cell = mesh.cell(ci);
        a[ai] = 1.0 / M_PI * exp(-(magSqr(cell.C() - Vector(0.25,0,0))));
    }

    auto grad = fvm::grad(a);

    a.output("a");
    auto S(a);
    forAll(S, si) {
        const auto& ci = a.index()[si];
        const auto& cell = mesh.cell(ci);
        S[si] = - grad[si] * Vector(1, 0, 0);
    }
    S.output("source");

    MPI::Finalize();

    /// run with no params
    logger.warn << "type \"./meso-mpi -help\" to get help." << std::endl;
    return 0;
}
