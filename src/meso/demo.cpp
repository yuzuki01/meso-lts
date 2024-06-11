#include "meso.h"


int demo_omp(int *p_argc, char ***p_argv, MESO::ArgParser &parser) {
    using namespace MESO;

    MPI::Initialize(p_argc, p_argv);

    auto mesh = Mesh::load_gambit("./cavity.neu");
    auto dvs_mesh = Solver::generate_newton_cotes(mesh.dimension(), 3, 50, 4.0);
    mesh.info();
    auto mpi_task = MPI::DVS_partition(dvs_mesh);
    dvs_mesh.info();

    auto f_cell = Solver::DistributionFunction(dvs_mesh.NCELL, Field<Scalar>(mesh, cell_field_flag));
    Field<Scalar> rho_cell(mesh, cell_field_flag);
    Field<Vector> vel_cell(mesh, cell_field_flag);
    auto f_maxwell = [](double density, double temperature, double cc) {
        return density / (M_PI * temperature) * exp(-cc / temperature);
    };

    double Rho0 = 1.0, T0 = 1.0;
    Vector U0(0.0, 0.0, 0.0);

    Field<Scalar> m0_local(mesh, cell_field_flag), m0_global(mesh, cell_field_flag);
    Field<Vector> m1_local(mesh, cell_field_flag), m1_global(mesh, cell_field_flag);

    clock_t s, e;

    s = clock();
    for (auto &cell: mesh.cells) {
        double m0 = 0.0;
        Vector m1(0.0, 0.0, 0.0);
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            auto c = particle.position - U0;
            auto cc = c * c;
            auto f = f_maxwell(Rho0, T0, cc);
            f_cell[p][cell.id] = f;
            m0 += particle.volume * f;
            m1 += particle.volume * f * particle.position;
        }
        m0_local[cell.id] = m0;
        m1_local[cell.id] = m1;
    }
    e = clock();
    logger.note << "Cost: " << double(e - s) / CLOCKS_PER_SEC << std::endl;

    MPI::AllReduce(m0_local, m0_global);
    MPI::AllReduce(m1_local, m1_global);

    for (auto &cell : mesh.cells) {
        auto rho = m0_global[cell.id];
        auto rhoU = m1_global[cell.id];
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = rhoU / rho;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    auto Uf = vel_cell.heft(0);
    auto Vf = vel_cell.heft(1);
    mesh.output("demo-omp",
                {"rho", "U", "V"},
                {&rho_cell, &Uf, &Vf});
    return 0;
}
