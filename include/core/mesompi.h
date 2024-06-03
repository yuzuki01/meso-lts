#ifndef CORE_MESOMPI_H
#define CORE_MESOMPI_H


namespace MESO::Math {
    /// defined in core/math.h
    class Vector;
}

namespace MESO {
    /// defined in mesh/object.h
    template<class FieldType>
    class Field;
}

namespace MESO::MPI {

    const int main_rank = 0;
    extern int node_num;
    extern int rank;
    extern int omp_num;

    struct MPI_TaskObject {
        int start = 0;
        int size = 0;
    };
    typedef std::vector<MPI_TaskObject> MPI_Task;
    MPI_Task get_task_distribution(int total_num);

    void Initialize(int *p_argc, char*** p_argv);

    void Finalize();

    void Bcast(Scalar &global);
    void Bcast(MESO::Vector &global);
    void ReduceAll(Scalar local, Scalar &global);
    void ReduceAll(const MESO::Vector &local, MESO::Vector &global);

    /// defined in solver/field.cpp
    void ReduceAll(Field<Scalar> &local, Field<Scalar> &global);
    void ReduceAll(Field<Vector> &local, Field<Vector> &global);
}

#endif //CORE_MESOMPI_H
