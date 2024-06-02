#ifndef CORE_MESOMPI_H
#define CORE_MESOMPI_H


namespace MESO::Math {
    /// defined in core/math.h
    class Vector;
}

namespace MESO {
    struct MPI_Task {
        size_t start = 0;
        size_t size = 0;
    };
    typedef std::vector<MPI_Task> MPI_Task_List;

    /// defined in mesh/object.h
    template<class FieldType>
    class Field;
}

namespace MESO::MPI {

    const int main_rank = 0;
    extern int processor_num;
    extern int rank;

    void Initialize(int *p_argc, char*** p_argv);

    void Finalize();

    void ReduceAndBcast(Scalar local, Scalar &global);
    void ReduceAndBcast(const MESO::Vector &local, MESO::Vector &global);
}

#endif //CORE_MESOMPI_H
