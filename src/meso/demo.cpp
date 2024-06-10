#include "meso.h"


int demo_omp(int *p_argc, char ***p_argv, MESO::ArgParser &parser) {
    using namespace MESO;
    auto mesh = Mesh::load_gambit("./cavity.neu");
    auto dvs_mesh = Solver::generate_newton_cotes(mesh.dimension(), 3, 33, 4.0);
    mesh.info();
    dvs_mesh.info();


    MPI::Initialize(p_argc, p_argv, parser.parse_param<int>("parallel", 1, true));

    auto mpi_task = MPI::DVS_partition(dvs_mesh);
    auto f_cell = Solver::DistributionFunction(dvs_mesh.NCELL, Field<Scalar>(mesh, cell_field_flag));
    Field<Scalar> rho_cell(mesh, cell_field_flag);
    Field<Vector> vel_cell(mesh, cell_field_flag);
    auto f_maxwell = [](double density, double temperature, double cc) {
        return density / (M_PI * temperature) * exp(-cc / temperature);
    };

    double Rho0 = 1.0, T0 = 1.0;
    Vector U(0.0, 0.0, 0.0);

    clock_t start = clock();

    for (auto &cell: mesh.cells) {
        double m0_local = 0.0;
        Vector m1_local(0.0, 0.0, 0.0);
        clock_t omp_start = clock();
#pragma omp parallel for shared(f_cell, cell, mpi_task, dvs_mesh, Rho0, T0, U, f_maxwell) \
        reduction(+:m0_local) reduction(+:m1_local) default(none)
        for (int p = 0; p < mpi_task.size; ++p) {
            ObjectId dvs_id = p + mpi_task.start;
            auto &particle = dvs_mesh.cells[dvs_id];
            auto c = particle.position - U;
            auto cc = c * c;
            auto f = f_maxwell(Rho0, T0, cc);
            f_cell[p][cell.id] = f;
            m0_local += particle.volume * f;
            m1_local += particle.volume * f * particle.position;
        }
        clock_t omp_end = clock();
#pragma omp critical
        {
            std::cout << " cell-" << cell.id << " omp cost: " << std::setprecision(8) << (double(omp_end - omp_start) / double(CLOCKS_PER_SEC)) << std::endl;
        }
        double m0;
        Vector m1;
        MPI::AllReduce(m0_local, m0);
        MPI::AllReduce(m1_local, m1);
        auto u = m1 / m0;
        rho_cell[cell.id] = m0;
        vel_cell[cell.id] = u;
    }
    clock_t end = clock();
    logger.note << "Cost: " << double(end - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    auto Uf = vel_cell.heft(0);
    auto Vf = vel_cell.heft(1);
    mesh.output("demo-omp",
                {"rho", "U", "V"},
                {&rho_cell, &Uf, &Vf});
    return 0;
}
